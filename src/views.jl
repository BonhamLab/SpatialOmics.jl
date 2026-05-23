# ── SpatialExtent ─────────────────────────────────────────────────────────────

"""
    SpatialExtent(xmin, xmax, ymin, ymax; coord_system="")
    SpatialExtent(shp::SpatialShapes)

Axis-aligned bounding box in a named coordinate system.

The `coord_system` field is checked when combining two extents (`union`,
`intersect`) or filtering elements with `view` — a mismatch raises an error
rather than silently producing wrong results.

```julia
ext = SpatialExtent(1000.0, 2000.0, 500.0, 1500.0; coord_system="global")
roi = view(ds, ext)
```

# See also
[`SpatialROI`](@ref), [`SpatialDatasetView`](@ref), [`SpatialElementView`](@ref)
"""
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

function SpatialExtent(shp::SpatialShapes)
    xmin = ymin =  Inf
    xmax = ymax = -Inf
    for g in shp.geometries
        for ring in GeoInterface.coordinates(g), pt in ring
            xmin = min(xmin, Float64(pt[1])); xmax = max(xmax, Float64(pt[1]))
            ymin = min(ymin, Float64(pt[2])); ymax = max(ymax, Float64(pt[2]))
        end
    end
    SpatialExtent(xmin, xmax, ymin, ymax; coord_system=shp.coord_system)
end

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

"""
    SpatialROI(geometry; coord_system="")

Polygon region of interest backed by any GeoInterface-compatible geometry.

Wraps the geometry with a cached axis-aligned bounding box for fast pre-filtering.
Point containment uses `GeometryOps.contains`; `pt ∈ roi` and
`filter(pt -> pt ∈ roi, pts)` are the idiomatic containment APIs.

# See also
[`SpatialExtent`](@ref), [`geometry`](@ref), [`SpatialElementView`](@ref)
"""
struct SpatialROI{G}
    geometry     :: G
    extent       :: SpatialExtent    # cached bbox of the geometry
    coord_system :: String
end

function SpatialROI(geom; coord_system::String="")
    SpatialROI(geom, _extent_of(geom, coord_system), coord_system)
end

"""
    geometry(roi) → geometry

Return the underlying GeoInterface geometry from a `SpatialROI`.
"""
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

# ── ROI → SpatialShapes ───────────────────────────────────────────────────────
# Symmetric with SpatialExtent(shp::SpatialShapes).

function SpatialShapes(ext::SpatialExtent; instance_id::Int32=Int32(1))
    ring = [Point2f(ext.xmin, ext.ymin), Point2f(ext.xmax, ext.ymin),
            Point2f(ext.xmax, ext.ymax), Point2f(ext.xmin, ext.ymax),
            Point2f(ext.xmin, ext.ymin)]
    SpatialShapes([Polygon(ring)]; instance_id=[instance_id], coord_system=ext.coord_system)
end

function SpatialShapes(roi::SpatialROI; instance_id::Int32=Int32(1))
    ring = Vector{Point2f}(GeoInterface.coordinates(roi.geometry)[1])
    first(ring) ≈ last(ring) || push!(ring, ring[1])
    SpatialShapes([Polygon(ring)]; instance_id=[instance_id], coord_system=roi.coord_system)
end

"""
    select(ax, ds) → SpatialROI

Interactively draw a polygon ROI on Makie axis `ax` and return the result as a
`SpatialROI`. Requires a Makie backend to be loaded.

This function is a stub; the implementation is provided by the Makie extension
when a backend is loaded.
"""
function select end

# ── SpatialElementView ────────────────────────────────────────────────────────

"""
    SpatialElementView{T, R}

Lazy spatial view into a single element, analogous to Julia's `SubArray`.

Holds a reference to the parent element and a region (`SpatialExtent` or
`SpatialROI`). The spatial filter is applied only when accessors (`coords`,
`geometries`, `feature_ids`, etc.) or `collect` are called — no data is copied
on construction.

`overlap` controls how shapes are matched:
- `:any` (default) — include shapes whose bounding box intersects the ROI
- `:full` — include only shapes fully contained within the ROI

# See also
[`SpatialDatasetView`](@ref), [`SpatialExtent`](@ref), [`SpatialROI`](@ref)
"""
struct SpatialElementView{T, R}
    parent  :: T
    roi     :: R
    overlap :: Symbol   # :any — shape intersects ROI; :full — shape fully inside ROI
end

# ── SpatialDatasetView ────────────────────────────────────────────────────────

"""
    SpatialDatasetView

Lazy view across all elements of a `SpatialDataset`, scoped to a spatial region.

Produced by `view(ds, extent)` or `view(ds, roi)`. Accessing an element via
the typed accessors (`points`, `shapes`, `images`, `labels`) returns a
`SpatialElementView` for that element — still lazy, no data copied.

```julia
roi = view(ds, SpatialExtent(1000.0, 2000.0, 500.0, 1500.0))
tx  = points(roi, "transcripts")   # SpatialElementView{SpatialPoints}
collect(tx)                        # materialise into a concrete SpatialPoints
```

# See also
[`SpatialElementView`](@ref), [`SpatialExtent`](@ref)
"""
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

# Shapes × SpatialExtent: per-geometry extent check, exact for rectangular ROIs
function _mask(shp::SpatialShapes, ext::SpatialExtent, overlap::Symbol=:any)
    if overlap == :any
        BitVector(let e = _extent_of(shp.geometries[i])
                  e.xmin <= ext.xmax && e.xmax >= ext.xmin &&
                  e.ymin <= ext.ymax && e.ymax >= ext.ymin
                  end for i in eachindex(shp.geometries))
    else  # :full — shape extent must lie entirely within ROI extent
        BitVector(let e = _extent_of(shp.geometries[i])
                  e.xmin >= ext.xmin && e.xmax <= ext.xmax &&
                  e.ymin >= ext.ymin && e.ymax <= ext.ymax
                  end for i in eachindex(shp.geometries))
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

function coords(v::SpatialElementView{<:SpatialPoints}, feature::String)
    p   = v.parent
    idx = findfirst(==(feature), p.feature_codebook)
    idx === nothing && return Point{2, eltype(eltype(p.coords))}[]
    mask = _mask(p, v.roi, v.overlap)
    p.coords[mask .& (p.feature_id .== idx)]
end

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
