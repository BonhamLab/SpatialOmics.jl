
# ── SpatialPoints ─────────────────────────────────────────────────────────────

mutable struct SpatialPoints{T<:AbstractFloat}
    coords           :: Vector{Point{2,T}}
    feature_id       :: Vector{Int32}
    feature_codebook :: Vector{String}
    instance_id      :: Vector{Int32}
    coord_system     :: String
    _attachment      :: Union{Nothing, Tuple{WeakRef, String}}
end

# Bare coords constructor
function SpatialPoints(coords::Vector{Point{2,T}};
                       feature_id::Vector{Int32}=zeros(Int32, length(coords)),
                       feature_codebook::Vector{String}=String[],
                       instance_id::Vector{Int32}=zeros(Int32, length(coords)),
                       coord_system::String="") where T<:AbstractFloat
    SpatialPoints{T}(coords, feature_id, feature_codebook, instance_id, coord_system, nothing)
end

# Tables.jl constructor — columns must have x and y; gene is optional
function SpatialPoints(table;
                       x::Symbol=:x, y::Symbol=:y,
                       gene::Union{Symbol,Nothing}=nothing,
                       coord_system::String="")
    cols = Tables.columntable(table)
    xs = cols[x]
    ys = cols[y]
    n = length(xs)
    T = Float32
    coords = [Point{2,T}(xs[i], ys[i]) for i in 1:n]
    if gene !== nothing && hasproperty(cols, gene)
        genes = cols[gene]
        codebook = unique(String.(genes))
        gene_idx = Dict(g => Int32(i) for (i, g) in enumerate(codebook))
        feature_id = [gene_idx[String(g)] for g in genes]
    else
        codebook = String[]
        feature_id = zeros(Int32, n)
    end
    SpatialPoints{T}(coords, feature_id, codebook, zeros(Int32, n), coord_system, nothing)
end

Base.length(pts::SpatialPoints) = length(pts.coords)

# Public accessors — field names are implementation detail
coords(pts::SpatialPoints)       = pts.coords
features(pts::SpatialPoints)     = pts.feature_codebook
coord_system(pts::SpatialPoints) = pts.coord_system
feature_ids(pts::SpatialPoints)  = pts.feature_id
instance_id(pts::SpatialPoints)  = pts.instance_id

# Filtered coords by feature name — hides the integer-index encoding from users
function coords(pts::SpatialPoints, feature::String)
    idx = findfirst(==(feature), pts.feature_codebook)
    idx === nothing && error("Feature \"$feature\" not in codebook")
    pts.coords[pts.feature_id .== idx]
end

# ── GeoInterface — MultiPoint ─────────────────────────────────────────────────

GeoInterface.isgeometry(::Type{<:SpatialPoints}) = true
GeoInterface.geomtrait(::SpatialPoints) = GeoInterface.MultiPointTrait()
GeoInterface.ngeom(::GeoInterface.MultiPointTrait, pts::SpatialPoints) = length(pts.coords)
GeoInterface.getgeom(::GeoInterface.MultiPointTrait, pts::SpatialPoints, i::Int) = pts.coords[i]

# ── SpatialShapes ─────────────────────────────────────────────────────────────

mutable struct SpatialShapes{G<:AbstractGeometry}
    geometries   :: Vector{G}
    bbox         :: Matrix{Float64}    # N×4 [xmin xmax ymin ymax]
    instance_id  :: Vector{Int32}
    coord_system :: String
    _attachment  :: Union{Nothing, Tuple{WeakRef, String}}
end

function _compute_bbox(geoms::Vector{<:AbstractGeometry})
    n = length(geoms)
    bbox = zeros(Float64, n, 4)
    for (i, g) in enumerate(geoms)
        rings = GeoInterface.coordinates(g)    # [[ring1_coords...], ...]
        xmin = ymin =  Inf
        xmax = ymax = -Inf
        for ring in rings, pt in ring
            xmin = min(xmin, Float64(pt[1]))
            xmax = max(xmax, Float64(pt[1]))
            ymin = min(ymin, Float64(pt[2]))
            ymax = max(ymax, Float64(pt[2]))
        end
        bbox[i, 1] = xmin; bbox[i, 2] = xmax
        bbox[i, 3] = ymin; bbox[i, 4] = ymax
    end
    bbox
end

function SpatialShapes(geometries::Vector{G};
                       instance_id::Vector{Int32}=zeros(Int32, length(geometries)),
                       coord_system::String="") where G<:AbstractGeometry
    bbox = _compute_bbox(geometries)
    SpatialShapes{G}(geometries, bbox, instance_id, coord_system, nothing)
end

Base.length(shp::SpatialShapes) = length(shp.geometries)

# ── Row type ──────────────────────────────────────────────────────────────────

struct SpatialShape{G<:AbstractGeometry}
    geometry     :: G
    instance_id  :: Int32
    bbox         :: NTuple{4, Float64}   # (xmin, xmax, ymin, ymax)
    coord_system :: String
end

function Base.getindex(shp::SpatialShapes{G}, i::Int) where G
    SpatialShape{G}(shp.geometries[i],
                    shp.instance_id[i],
                    (shp.bbox[i,1], shp.bbox[i,2], shp.bbox[i,3], shp.bbox[i,4]),
                    shp.coord_system)
end

Base.iterate(shp::SpatialShapes, i=1) = i > length(shp) ? nothing : (shp[i], i+1)
Base.eltype(::Type{SpatialShapes{G}}) where G = SpatialShape{G}

function Base.filter(pred, shp::SpatialShapes{G}) where G
    keep = [i for i in eachindex(shp.geometries) if pred(shp[i])]
    SpatialShapes(shp.geometries[keep];
                  instance_id  = shp.instance_id[keep],
                  coord_system = shp.coord_system)
end

# Accessors
geometries(shp::SpatialShapes)  = shp.geometries
bbox(shp::SpatialShapes)        = shp.bbox
instance_id(shp::SpatialShapes) = shp.instance_id
instance_id(s::SpatialShape)    = s.instance_id
coord_system(shp::SpatialShapes) = shp.coord_system

# ── GeoInterface — GeometryCollection ─────────────────────────────────────────

GeoInterface.isgeometry(::Type{<:SpatialShapes}) = true
GeoInterface.geomtrait(::SpatialShapes) = GeoInterface.GeometryCollectionTrait()
GeoInterface.ngeom(::GeoInterface.GeometryCollectionTrait, shp::SpatialShapes) = length(shp.geometries)
GeoInterface.getgeom(::GeoInterface.GeometryCollectionTrait, shp::SpatialShapes, i::Int) = shp.geometries[i]

# ── Geometry transform helper — extend for other types as needed ───────────────

function _transform_geom(t::AbstractTransformation, poly::Polygon)
    rings = GeoInterface.coordinates(poly)    # [[exterior_pts...], [hole_pts...], ...]
    T = Float32
    new_rings = map(rings) do ring
        [let v = apply(t, SVector{2,Float64}(pt[1], pt[2])); Point{2,T}(v[1], v[2]); end
         for pt in ring]
    end
    ext = new_rings[1]
    holes = length(new_rings) > 1 ? new_rings[2:end] : Vector{Vector{Point{2,T}}}()
    Polygon(ext, holes)
end

# ── apply / apply! on SpatialPoints ──────────────────────────────────────────

function apply(t::AbstractTransformation, pts::SpatialPoints{T}) where T
    new_coords = map(pts.coords) do p
        v = apply(t, p)
        Point{2,T}(v[1], v[2])
    end
    SpatialPoints{T}(new_coords, copy(pts.feature_id), copy(pts.feature_codebook),
                     copy(pts.instance_id), t.dst, nothing)
end

function apply!(t::AbstractTransformation, pts::SpatialPoints{T}) where T
    map!(pts.coords, pts.coords) do p
        v = apply(t, p)
        Point{2,T}(v[1], v[2])
    end
    pts.coord_system = t.dst
    pts
end

# ── apply / apply! on SpatialShapes ───────────────────────────────────────────

function apply(t::AbstractTransformation, shp::SpatialShapes{G}) where G
    new_geoms = G[_transform_geom(t, g) for g in shp.geometries]
    SpatialShapes(new_geoms; instance_id=copy(shp.instance_id), coord_system=t.dst)
end

function apply!(t::AbstractTransformation, shp::SpatialShapes{G}) where G
    for i in eachindex(shp.geometries)
        shp.geometries[i] = _transform_geom(t, shp.geometries[i])
    end
    shp.bbox = _compute_bbox(shp.geometries)
    shp.coord_system = t.dst
    shp
end

# ── Back-reference: element ↔ dataset ownership ───────────────────────────────
#
# Elements carry a WeakRef to their dataset so dispatch can resolve mappings
# without the dataset being passed explicitly. WeakRef prevents elements from
# keeping a dataset alive past its scope.
#
# Invariant: each element belongs to at most one dataset at a time.
# setindex!(ds, el, name) enforces this; copy(el) strips the reference.

_dataset_ref(el::SpatialPoints)  = el._attachment
_dataset_ref(el::SpatialShapes)  = el._attachment
_dataset_ref(::Any)              = nothing   # images, labels, tables: no ref yet

function _owning_dataset(el)
    att = _dataset_ref(el)
    att === nothing && return nothing
    att[1].value   # WeakRef → live dataset or nothing if GC'd
end

function _set_backref!(el::Union{SpatialPoints, SpatialShapes},
                       ds::SpatialDataset, name::String)
    el._attachment = (WeakRef(ds), name)
end
_set_backref!(::Any, ::SpatialDataset, ::String) = nothing  # no-op for other types

function Base.copy(pts::SpatialPoints{T}) where T
    SpatialPoints{T}(copy(pts.coords), copy(pts.feature_id), copy(pts.feature_codebook),
                     copy(pts.instance_id), pts.coord_system, nothing)
end

function Base.copy(shp::SpatialShapes{G}) where G
    SpatialShapes{G}(copy(shp.geometries), copy(shp.bbox), copy(shp.instance_id),
                     shp.coord_system, nothing)
end

# ── Typed dataset accessors ───────────────────────────────────────────────────

function points(ds::SpatialDataset, name::String)
    el = ds.elements[name]
    el isa SpatialPoints || error("Element \"$name\" is not SpatialPoints (got $(typeof(el)))")
    el
end

function shapes(ds::SpatialDataset, name::String)
    el = ds.elements[name]
    el isa SpatialShapes || error("Element \"$name\" is not SpatialShapes (got $(typeof(el)))")
    el
end
