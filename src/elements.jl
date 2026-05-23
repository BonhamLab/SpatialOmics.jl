
# ── SpatialPoints ─────────────────────────────────────────────────────────────

"""
    SpatialPoints{T<:AbstractFloat}

Point cloud of 2-D spatial observations — transcripts, centroids, or other
coordinates in a named coordinate system.

Stored as parallel arrays (struct-of-arrays layout). Feature labels (gene names,
cell types, etc.) are encoded as integer indices into a compact `feature_codebook`
for O(1) lookup by name via `coords(pts, feature)`.

# Constructors

    SpatialPoints(coords; feature_id, feature_codebook, instance_id, coord_system)

Bare coordinates constructor. `coords` is a `Vector{Point{2,T}}`.

    SpatialPoints(table; x=:x, y=:y, gene=nothing, coord_system="")

Tables.jl constructor. Reads x/y from columns named by `x` and `y`; optionally
encodes a gene/label column via `gene`.

```julia
pts = SpatialPoints(df; x=:x_centroid, y=:y_centroid, gene=:target, coord_system="global")
length(pts)          # number of points
features(pts)        # gene names in codebook order
coords(pts, "Epcam") # coordinates of all Epcam transcripts
```

# See also
[`coords`](@ref), [`features`](@ref), [`feature_ids`](@ref),
[`instance_id`](@ref), [`subsample`](@ref), [`top_features`](@ref)
"""
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

"""
    coords(pts) → Vector{Point{2,T}}
    coords(pts, feature) → Vector{Point{2,T}}

Return all coordinates in `pts`, or only those whose feature label matches `feature`.

The single-argument form returns the full coordinate vector. The two-argument
form performs an O(n) filter using the integer codebook index; raises an error
if `feature` is not in the codebook.

# See also
[`features`](@ref), [`feature_ids`](@ref)
"""
coords(pts::SpatialPoints)       = pts.coords

"""
    features(pts) → Vector{String}

Return the feature codebook — the unique feature labels (gene names, cell types, …)
present in `pts`, in the order used by `feature_ids`.
"""
features(pts::SpatialPoints)     = pts.feature_codebook

"""
    coord_system(el) → String

Return the name of the coordinate system that `el` belongs to.

Defined for `SpatialPoints`, `SpatialShapes`, `SpatialImage`, `SpatialLabels`,
`SpatialExtent`, `SpatialROI`, and their view types.
"""
coord_system(pts::SpatialPoints) = pts.coord_system

"""
    feature_ids(pts) → Vector{Int32}

Return the per-point feature index vector. Each value is a 1-based index into
`features(pts)`; `0` means unlabelled.

# See also
[`features`](@ref), [`coords`](@ref)
"""
feature_ids(pts::SpatialPoints)  = pts.feature_id

"""
    instance_id(el) → Vector{Int32}
    instance_id(shape) → Int32

Return per-observation instance IDs for a collection, or the single instance
ID for a `SpatialShape` row. `0` means unassigned (transcript not inside any
cell, or shape not assigned to a region).

# See also
[`instance_ids`](@ref), [`count_per_instance`](@ref)
"""
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

"""
    SpatialShapes{G<:AbstractGeometry}

Collection of 2-D geometries (typically cell boundary polygons) in a named
coordinate system. Each shape carries an `instance_id` linking it to cells
or objects in a `SpatialRelation`.

Implements the GeoInterface `GeometryCollectionTrait`, making it compatible
with GeometryOps operations directly.

# Constructors
    SpatialShapes(geometries; instance_id, coord_system)

    SpatialShapes(ext::SpatialExtent)        # rectangular region
    SpatialShapes(roi::SpatialROI)           # polygon region

# See also
[`geometries`](@ref), [`instance_id`](@ref), [`SpatialROI`](@ref), [`SpatialPoints`](@ref)
"""
mutable struct SpatialShapes{G<:AbstractGeometry}
    geometries   :: Vector{G}
    instance_id  :: Vector{Int32}
    coord_system :: String
    _attachment  :: Union{Nothing, Tuple{WeakRef, String}}
end

function SpatialShapes(geometries::Vector{G};
                       instance_id::Vector{Int32}=zeros(Int32, length(geometries)),
                       coord_system::String="") where G<:AbstractGeometry
    SpatialShapes{G}(geometries, instance_id, coord_system, nothing)
end

Base.length(shp::SpatialShapes) = length(shp.geometries)

# ── Row type ──────────────────────────────────────────────────────────────────

"""
    SpatialShape{G<:AbstractGeometry}

Single-shape row accessor produced by indexing into a `SpatialShapes` collection.

Carries the geometry, its `instance_id`, and the coordinate system name.
Row-accessor and collection share the same field names (`geometry`, `instance_id`,
`coord_system`) so code generalises across both.

# See also
[`SpatialShapes`](@ref), [`geometry`](@ref)
"""
struct SpatialShape{G<:AbstractGeometry}
    geometry     :: G
    instance_id  :: Int32
    coord_system :: String
end

Base.getindex(shp::SpatialShapes{G}, i::Int) where G =
    SpatialShape{G}(shp.geometries[i], shp.instance_id[i], shp.coord_system)

Base.iterate(shp::SpatialShapes, i=1) = i > length(shp) ? nothing : (shp[i], i+1)
Base.eltype(::Type{SpatialShapes{G}}) where G = SpatialShape{G}

function Base.filter(pred, shp::SpatialShapes{G}) where G
    keep = [i for i in eachindex(shp.geometries) if pred(shp[i])]
    SpatialShapes(shp.geometries[keep];
                  instance_id  = shp.instance_id[keep],
                  coord_system = shp.coord_system)
end

"""
    geometries(shp) → Vector{G}

Return the vector of geometries from a `SpatialShapes` collection.

# See also
[`SpatialShapes`](@ref), [`geometry`](@ref)
"""
geometries(shp::SpatialShapes)   = shp.geometries
instance_id(shp::SpatialShapes)  = shp.instance_id
instance_id(s::SpatialShape)     = s.instance_id
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

"""
    apply!(t, el) → el

Apply transformation `t` to `el` in-place, mutating coordinates and updating
the element's `coord_system` to `t.dst`. Returns `el`.

Defined for `SpatialPoints` and `SpatialShapes`. Prefer `apply` (non-mutating)
when the element is attached to a dataset.

# See also
[`apply`](@ref)
"""
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

Base.copy(shp::SpatialShapes{G}) where G =
    SpatialShapes{G}(copy(shp.geometries), copy(shp.instance_id), shp.coord_system, nothing)

"""
    instance_ids(pts) → Vector{Int32}

Return the per-point instance ID vector for a `SpatialPoints` collection.

Alias for `instance_id(pts)` provided for consistency with
`instance_ids(lbl::SpatialLabels)`. `0` means unassigned.
"""
instance_ids(pts::SpatialPoints) = pts.instance_id

# ── subsample ─────────────────────────────────────────────────────────────────

"""
    subsample(pts, n) → SpatialPoints

Return a new `SpatialPoints` with a random subset of `n` observations.

If `n ≥ length(pts)`, the original object is returned unchanged.
The feature codebook is preserved; indices are rebuilt from the subset.

# See also
[`top_features`](@ref)
"""
function subsample(pts::SpatialPoints{T}, n::Int) where T
    n >= length(pts) && return pts
    idx = sort!(randperm(length(pts))[1:n])
    SpatialPoints{T}(pts.coords[idx], pts.feature_id[idx], copy(pts.feature_codebook),
                     pts.instance_id[idx], pts.coord_system, nothing)
end

# ── top_features ──────────────────────────────────────────────────────────────

"""
    top_features(pts, n=10) → Vector{String}

Return the `n` most frequent feature names in `pts`, sorted by descending count.

Returns an empty vector if `pts` has no feature codebook.

# See also
[`count_per_instance`](@ref), [`features`](@ref)
"""
function top_features(pts::SpatialPoints, n::Int=10)
    cb = pts.feature_codebook
    isempty(cb) && return String[]
    counts = [count(==(Int32(i)), pts.feature_id) for i in eachindex(cb)]
    idx    = sortperm(counts; rev=true)[1:min(n, length(cb))]
    cb[idx]
end

# ── count_per_instance ────────────────────────────────────────────────────────

"""
    count_per_instance(pts) → Dict{Int32, Int}

Return a dictionary mapping each non-zero instance ID to its observation count.

Unassigned points (`instance_id == 0`) are excluded. Useful for computing
transcript counts per cell or density metrics.

# See also
[`instance_id`](@ref), [`top_features`](@ref)
"""
function count_per_instance(pts::SpatialPoints)
    counts = Dict{Int32, Int}()
    for id in pts.instance_id
        id == Int32(0) && continue
        counts[id] = get(counts, id, 0) + 1
    end
    counts
end

# ── Typed dataset accessors ───────────────────────────────────────────────────

"""
    points(ds, name) → SpatialPoints
    points(v, name) → SpatialElementView{SpatialPoints}

Retrieve the named `SpatialPoints` element from a dataset or dataset view.
Raises an error if the element exists but is not `SpatialPoints`.

# See also
[`shapes`](@ref), [`images`](@ref), [`labels`](@ref), [`elements`](@ref)
"""
function points(ds::SpatialDataset, name::String)
    el = ds.elements[name]
    el isa SpatialPoints || error("Element \"$name\" is not SpatialPoints (got $(typeof(el)))")
    el
end

"""
    shapes(ds, name) → SpatialShapes
    shapes(v, name) → SpatialElementView{SpatialShapes}

Retrieve the named `SpatialShapes` element from a dataset or dataset view.
Raises an error if the element exists but is not `SpatialShapes`.

# See also
[`points`](@ref), [`images`](@ref), [`labels`](@ref)
"""
function shapes(ds::SpatialDataset, name::String)
    el = ds.elements[name]
    el isa SpatialShapes || error("Element \"$name\" is not SpatialShapes (got $(typeof(el)))")
    el
end
