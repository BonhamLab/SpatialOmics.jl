# ── SpatialExtent ─────────────────────────────────────────────────────────────

struct SpatialExtent
    xmin :: Float64
    xmax :: Float64
    ymin :: Float64
    ymax :: Float64
    coord_system :: String
end

SpatialExtent(xmin, xmax, ymin, ymax; coord_system::String="") =
    SpatialExtent(Float64(xmin), Float64(xmax), Float64(ymin), Float64(ymax), coord_system)

coord_system(ext::SpatialExtent) = ext.coord_system

SpatialExtent(shp::SpatialShapes) =
    SpatialExtent(minimum(shp.bbox[:,1]), maximum(shp.bbox[:,2]),
                  minimum(shp.bbox[:,3]), maximum(shp.bbox[:,4]);
                  coord_system = shp.coord_system)

function _check_cs(a::SpatialExtent, b::SpatialExtent)
    a.coord_system == b.coord_system ||
        error("coord_system mismatch: \"$(a.coord_system)\" vs \"$(b.coord_system)\"")
end

function Base.union(a::SpatialExtent, b::SpatialExtent)
    _check_cs(a, b)
    SpatialExtent(min(a.xmin, b.xmin), max(a.xmax, b.xmax),
                  min(a.ymin, b.ymin), max(a.ymax, b.ymax);
                  coord_system = a.coord_system)
end

function Base.intersect(a::SpatialExtent, b::SpatialExtent)
    _check_cs(a, b)
    xmin, xmax = max(a.xmin, b.xmin), min(a.xmax, b.xmax)
    ymin, ymax = max(a.ymin, b.ymin), min(a.ymax, b.ymax)
    xmin < xmax && ymin < ymax || return nothing
    SpatialExtent(xmin, xmax, ymin, ymax; coord_system = a.coord_system)
end

# ── SpatialROI ────────────────────────────────────────────────────────────────

struct SpatialROI{G}
    geometry     :: G
    extent       :: SpatialExtent    # cached bbox of the geometry
    coord_system :: String
end

function SpatialROI(geom; coord_system::String="")
    SpatialROI(geom, _extent_of(geom, coord_system), coord_system)
end

geometry(roi::SpatialROI)      = roi.geometry
coord_system(roi::SpatialROI)  = roi.coord_system

function _extent_of(geom, cs::String="")
    rings = GeoInterface.coordinates(geom)
    xmin = ymin =  Inf
    xmax = ymax = -Inf
    for ring in rings, pt in ring
        xmin = min(xmin, Float64(pt[1])); xmax = max(xmax, Float64(pt[1]))
        ymin = min(ymin, Float64(pt[2])); ymax = max(ymax, Float64(pt[2]))
    end
    SpatialExtent(xmin, xmax, ymin, ymax; coord_system=cs)
end

# ── SpatialElementView ────────────────────────────────────────────────────────

struct SpatialElementView{T, R}
    parent  :: T
    roi     :: R
    overlap :: Symbol   # :any — shape intersects ROI; :full — shape fully inside ROI
end

# ── SpatialDatasetView ────────────────────────────────────────────────────────

struct SpatialDatasetView
    parent :: SpatialDataset
    roi    :: Union{SpatialExtent, SpatialROI}
end

# ── view constructors ─────────────────────────────────────────────────────────

const _ROI = Union{SpatialExtent, SpatialROI}

function Base.view(el::Union{SpatialPoints, SpatialShapes}, roi::_ROI;
                   overlap::Symbol=:any)
    overlap in (:any, :full) || error("overlap must be :any or :full, got :$overlap")
    cs_el  = coord_system(el)
    cs_roi = coord_system(roi)
    cs_el == cs_roi || error(
        "Coordinate system mismatch: element is \"$cs_el\", ROI is \"$cs_roi\". " *
        "Apply a transform to the element first.")
    SpatialElementView(el, roi, overlap)
end

function Base.view(ds::SpatialDataset, roi::_ROI)
    SpatialDatasetView(ds, roi)
end

# ── SpatialDatasetView element access ─────────────────────────────────────────

Base.getindex(v::SpatialDatasetView, name::String) = view(v.parent[name], v.roi)
Base.haskey(v::SpatialDatasetView, name::String)   = haskey(v.parent, name)
Base.keys(v::SpatialDatasetView)                   = keys(v.parent)

function points(v::SpatialDatasetView, name::String)
    el = v.parent.elements[name]
    el isa SpatialPoints || error("Element \"$name\" is not SpatialPoints (got $(typeof(el)))")
    view(el, v.roi)
end

function shapes(v::SpatialDatasetView, name::String)
    el = v.parent.elements[name]
    el isa SpatialShapes || error("Element \"$name\" is not SpatialShapes (got $(typeof(el)))")
    view(el, v.roi)
end

# ── Mask computation ──────────────────────────────────────────────────────────

# Points: overlap mode is irrelevant — a point is either inside or not
function _mask(pts::SpatialPoints, ext::SpatialExtent, ::Symbol=:any)
    BitVector(ext.xmin <= p[1] <= ext.xmax && ext.ymin <= p[2] <= ext.ymax
              for p in pts.coords)
end

function _mask(pts::SpatialPoints, roi::SpatialROI, ::Symbol=:any)
    pre = _mask(pts, roi.extent)
    mask = copy(pre)
    for i in findall(pre)
        mask[i] = GeometryOps.contains(roi.geometry, pts.coords[i])
    end
    mask
end

# Shapes × SpatialExtent: bbox-based check is exact for rectangular ROIs
function _mask(shp::SpatialShapes, ext::SpatialExtent, overlap::Symbol=:any)
    bb = shp.bbox   # N×4 [xmin xmax ymin ymax]
    if overlap == :any
        BitVector(bb[i,1] <= ext.xmax && bb[i,2] >= ext.xmin &&
                  bb[i,3] <= ext.ymax && bb[i,4] >= ext.ymin
                  for i in 1:size(bb, 1))
    else  # :full — shape bbox must lie entirely within extent
        BitVector(bb[i,1] >= ext.xmin && bb[i,2] <= ext.xmax &&
                  bb[i,3] >= ext.ymin && bb[i,4] <= ext.ymax
                  for i in 1:size(bb, 1))
    end
end

# Shapes × SpatialROI: bbox pre-filter always :any, then exact polygon test
function _mask(shp::SpatialShapes, roi::SpatialROI, overlap::Symbol=:any)
    pre = _mask(shp, roi.extent)   # :any pre-filter — never excludes candidates
    mask = copy(pre)
    for i in findall(pre)
        mask[i] = overlap == :any ?
            GeometryOps.intersects(roi.geometry, shp.geometries[i]) :
            GeometryOps.contains(roi.geometry, shp.geometries[i])
    end
    mask
end

# ── collect — materialises a lazy view into a new element ─────────────────────

function Base.collect(v::SpatialElementView{<:SpatialPoints})
    mask = _mask(v.parent, v.roi, v.overlap)
    p = v.parent
    T = eltype(eltype(p.coords))
    SpatialPoints{T}(p.coords[mask], p.feature_id[mask], copy(p.feature_codebook),
                     p.instance_id[mask], p.coord_system, nothing)
end

function Base.collect(v::SpatialElementView{<:SpatialShapes})
    mask = _mask(v.parent, v.roi, v.overlap)
    s = v.parent
    SpatialShapes(s.geometries[mask];
                  instance_id=s.instance_id[mask], coord_system=s.coord_system)
end

# ── length — count without allocating a copy ──────────────────────────────────

Base.length(v::SpatialElementView) = count(_mask(v.parent, v.roi, v.overlap))

# ── passthrough accessors for SpatialElementView ──────────────────────────────

coord_system(v::SpatialElementView) = coord_system(v.parent)
coord_system(v::SpatialDatasetView) = coord_system(v.roi)

features(v::SpatialElementView{<:SpatialPoints})     = v.parent.feature_codebook

geometries(v::SpatialElementView{<:SpatialShapes}) =
    v.parent.geometries[_mask(v.parent, v.roi, v.overlap)]

coords(v::SpatialElementView{<:SpatialPoints}) =
    v.parent.coords[_mask(v.parent, v.roi, v.overlap)]

feature_ids(v::SpatialElementView{<:SpatialPoints}) =
    v.parent.feature_id[_mask(v.parent, v.roi, v.overlap)]

instance_id(v::SpatialElementView{<:SpatialShapes}) =
    v.parent.instance_id[_mask(v.parent, v.roi, v.overlap)]

instance_id(v::SpatialElementView{<:SpatialPoints}) =
    v.parent.instance_id[_mask(v.parent, v.roi, v.overlap)]

function count_per_instance(v::SpatialElementView{<:SpatialPoints})
    mask   = _mask(v.parent, v.roi, v.overlap)
    counts = Dict{Int32, Int}()
    for (i, id) in enumerate(v.parent.instance_id)
        mask[i] || continue
        id == Int32(0) && continue
        counts[id] = get(counts, id, 0) + 1
    end
    counts
end

# ── typed accessors on SpatialDatasetView ─────────────────────────────────────

function images(v::SpatialDatasetView, name::String)
    el = v.parent.elements[name]
    el isa SpatialImage || error("Element \"$name\" is not SpatialImage (got $(typeof(el)))")
    ext = v.roi isa SpatialExtent ? v.roi : v.roi.extent
    Base.view(el, ext)
end

function labels(v::SpatialDatasetView, name::String)
    el = v.parent.elements[name]
    el isa SpatialLabels || error("Element \"$name\" is not SpatialLabels (got $(typeof(el)))")
    el    # labels are rasters — no spatial element view; return as-is
end

function tables(v::SpatialDatasetView, name::String)
    el = v.parent.elements[name]
    el isa SpatialTable || error("Element \"$name\" is not SpatialTable (got $(typeof(el)))")
    el    # table filtering is driven by feature(); return as-is
end

# ── feature API — ROI-filtered (values in shape order within ROI) ──────────────

function feature(dsv::SpatialDatasetView, gene::String; region::String)
    shape_el = dsv.parent.elements[region]
    shape_el isa SpatialShapes || error("Element \"$region\" is not SpatialShapes")
    mask        = _mask(shape_el, dsv.roi, :any)
    ids_ordered = shape_el.instance_id[mask]

    for (_, el) in dsv.parent.elements
        el isa SpatialTable || continue
        el.region == region || continue
        idx = findfirst(==(gene), var_names(el))
        idx === nothing && error("Gene \"$gene\" not found in table linked to \"$region\"")
        obs_ids   = _obs_instance_ids(el)
        id_to_row = Dict(id => i for (i, id) in enumerate(obs_ids))
        return [el.X[id_to_row[id], idx] for id in ids_ordered]
    end
    error("No SpatialTable linked to element \"$region\"")
end
