
# ── SpatialPoints ─────────────────────────────────────────────────────────────

function _origin_fields(n::Int, ::Nothing, ::Nothing, codebook)
    isempty(codebook) || throw(ArgumentError("origin_codebook requires origin_id"))
    nothing, String[]
end

function _origin_fields(n::Int, labels::AbstractVector, ::Nothing, codebook)
    isempty(codebook) || throw(ArgumentError(
        "origin_codebook cannot be combined with origins",
    ))
    length(labels) == n || throw(DimensionMismatch(
        "origins has length $(length(labels)); expected $n",
    ))
    names = unique(String.(labels))
    positions = Dict(name => Int32(i) for (i, name) in enumerate(names))
    Int32[positions[String(label)] for label in labels], names
end

function _origin_fields(n::Int, ::Nothing, ids::AbstractVector{<:Integer}, codebook)
    length(ids) == n || throw(DimensionMismatch(
        "origin_id has length $(length(ids)); expected $n",
    ))
    names = String.(codebook)
    encoded = Int32.(ids)
    all(id -> 1 <= id <= length(names), encoded) || throw(ArgumentError(
        "origin_id values must index origin_codebook",
    ))
    encoded, names
end

function _origin_fields(::Int, ::AbstractVector, ::AbstractVector, _)
    throw(ArgumentError("origins and origin_id are alternative inputs"))
end

_subset_origin_ids(::Nothing, _) = nothing
_subset_origin_ids(ids::Vector{Int32}, idx) = ids[idx]

"""
    SpatialPoints{T<:AbstractFloat}

Point cloud of 2-D spatial observations — transcripts, centroids, or other
coordinates in a named coordinate system.

Stored as parallel arrays (struct-of-arrays layout). Feature labels (gene names,
cell types, etc.) are encoded as integer indices into a compact `feature_codebook`
for O(1) lookup by name via `coords(pts, feature)`.

# Constructors

    SpatialPoints(coords; feature_id, feature_codebook, instance_id, origins, coord_system)

Bare coordinates constructor. `coords` is a `Vector{Point{2,T}}`. Pass source
names with `origins`, or an encoded `origin_id` vector and `origin_codebook`,
to retain acquisition provenance independently of position.

    SpatialPoints(table; x=:x, y=:y, gene=nothing, origin=nothing, coord_system="")

Tables.jl constructor. Reads x/y from columns named by `x` and `y`; optionally
encodes gene/label and acquisition-source columns via `gene` and `origin`.

```julia
pts = SpatialPoints(df; x=:x_centroid, y=:y_centroid, gene=:target, coord_system="global")
length(pts)          # number of points
features(pts)        # gene names in codebook order
coords(pts, "Epcam") # coordinates of all Epcam transcripts
```

# See also
[`coords`](@ref), [`features`](@ref), [`feature_ids`](@ref), [`origins`](@ref),
[`instance_id`](@ref), [`subsample`](@ref), [`top_features`](@ref)
"""
mutable struct SpatialPoints{T<:AbstractFloat}
    coords           :: Vector{Point{2,T}}
    feature_id       :: Vector{Int32}
    feature_codebook :: Vector{String}
    instance_id      :: Vector{Int32}
    feature_columns  :: Union{Nothing, NamedTuple}
    origin_id        :: Union{Nothing, Vector{Int32}}
    origin_codebook  :: Vector{String}
    coord_system     :: String
    _attachment      :: Union{Nothing, Tuple{WeakRef, String}}
end

# Bare coords constructor
function SpatialPoints(coords::Vector{Point{2,T}};
                       feature_id::Vector{Int32}=zeros(Int32, length(coords)),
                       feature_codebook::Vector{String}=String[],
                       instance_id::Vector{Int32}=zeros(Int32, length(coords)),
                       features::Union{Nothing, NamedTuple}=nothing,
                       origins::Union{Nothing,AbstractVector}=nothing,
                       origin_id::Union{Nothing,AbstractVector{<:Integer}}=nothing,
                       origin_codebook::AbstractVector{<:AbstractString}=String[],
                       coord_system::String="") where T<:AbstractFloat
    encoded_origins, origin_names = _origin_fields(
        length(coords), origins, origin_id, origin_codebook,
    )
    SpatialPoints{T}(
        coords, feature_id, feature_codebook, instance_id, features,
        encoded_origins, origin_names, coord_system, nothing,
    )
end

# Tables.jl constructor — columns must have x and y; gene and features are optional
function SpatialPoints(table;
                       x::Symbol=:x, y::Symbol=:y,
                       gene::Union{Symbol,Nothing}=nothing,
                       origin::Union{Symbol,Nothing}=nothing,
                       features::Union{Nothing, NamedTuple}=nothing,
                       coord_system::String="")
    cols = Tables.columntable(table)
    xs = cols[x]
    ys = cols[y]
    n = length(xs)
    coords = [Point{2,Float32}(xs[i], ys[i]) for i in 1:n]
    if gene !== nothing && hasproperty(cols, gene)
        genes = cols[gene]
        codebook = unique(String.(genes))
        gene_idx = Dict(g => Int32(i) for (i, g) in enumerate(codebook))
        feature_id = [gene_idx[String(g)] for g in genes]
    else
        codebook = String[]
        feature_id = zeros(Int32, n)
    end
    origin_values = origin === nothing ? nothing : cols[origin]
    SpatialPoints(
        coords;
        feature_id,
        feature_codebook=codebook,
        instance_id=zeros(Int32, n),
        features,
        origins=origin_values,
        coord_system,
    )
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
    features(pts, col) → Vector

Return the feature codebook (gene names / cell types) when called with one argument,
or a named extra column when called with a `Symbol`:

```julia
features(pts)     # → ["Actb", "Epcam", …]  gene names
features(pts, :z) # → [3, 5, 3, …]           z-slice per point
```

Extra columns are stored via the `features` keyword at construction time:
```julia
pts = SpatialPoints(df; gene=:gene, features=(z=df.z,), coord_system="global_px")
```
"""
features(pts::SpatialPoints)                = pts.feature_codebook
features(pts::SpatialPoints, col::Symbol)   = getproperty(pts.feature_columns, col)

_subset_feature_columns(::Nothing, _)       = nothing
_subset_feature_columns(nt::NamedTuple, idx) =
    NamedTuple{keys(nt)}(map(v -> v[idx], values(nt)))

"""
    with_instance_ids(points, ids) -> SpatialPoints

Return a detached copy of `points` with replacement instance assignments.
Coordinates, feature encodings, auxiliary feature columns, acquisition origins,
and the coordinate system are preserved. `ids` must contain one value per point.

This is useful when importing assignments from an external segmentation tool
without rebuilding a point collection field by field.
"""
function with_instance_ids(pts::SpatialPoints{T}, ids::AbstractVector{<:Integer}) where T
    length(ids) == length(pts) || throw(DimensionMismatch(
        "instance IDs have length $(length(ids)); expected $(length(pts))",
    ))
    SpatialPoints{T}(
        copy(pts.coords), copy(pts.feature_id), copy(pts.feature_codebook),
        Int32.(ids), _subset_feature_columns(pts.feature_columns, :),
        isnothing(pts.origin_id) ? nothing : copy(pts.origin_id),
        copy(pts.origin_codebook), pts.coord_system, nothing,
    )
end

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
    SpatialShapes(geometries; instance_id, origins, coord_system)

    SpatialShapes(ext::SpatialExtent)        # rectangular region
    SpatialShapes(roi::SpatialROI)           # polygon region

# See also
[`geometries`](@ref), [`instance_id`](@ref), [`origins`](@ref),
[`SpatialROI`](@ref), [`SpatialPoints`](@ref)
"""
mutable struct SpatialShapes{G<:AbstractGeometry}
    geometries   :: Vector{G}
    instance_id  :: Vector{Int32}
    origin_id    :: Union{Nothing, Vector{Int32}}
    origin_codebook :: Vector{String}
    coord_system :: String
    _attachment  :: Union{Nothing, Tuple{WeakRef, String}}
end

function SpatialShapes(geometries::Vector{G};
                       instance_id::Vector{Int32}=zeros(Int32, length(geometries)),
                       origins::Union{Nothing,AbstractVector}=nothing,
                       origin_id::Union{Nothing,AbstractVector{<:Integer}}=nothing,
                       origin_codebook::AbstractVector{<:AbstractString}=String[],
                       coord_system::String="") where G<:AbstractGeometry
    encoded_origins, origin_names = _origin_fields(
        length(geometries), origins, origin_id, origin_codebook,
    )
    SpatialShapes{G}(
        geometries, instance_id, encoded_origins, origin_names, coord_system, nothing,
    )
end

Base.length(shp::SpatialShapes) = length(shp.geometries)

"""
    origins(el) → Vector{String}

Return the acquisition-source codebook for a point or shape element.

Source provenance is independent of geometry. Use `view(ds, source_name)` to
select observations acquired by one source, including when source footprints
overlap.
"""
origins(el::Union{SpatialPoints,SpatialShapes}) = el.origin_codebook

"""
    origin_ids(el) → Union{Nothing,Vector{Int32}}

Return compact per-observation indices into [`origins`](@ref), or `nothing`
when the element has no acquisition provenance.
"""
origin_ids(el::Union{SpatialPoints,SpatialShapes}) = el.origin_id

"""
    source(el, i) → Union{Nothing,String}

Return the acquisition-source name for observation `i`, or `nothing` when the
element has no acquisition provenance.
"""
source(el::Union{SpatialPoints,SpatialShapes}, i::Integer) =
    el.origin_id === nothing ? nothing : el.origin_codebook[el.origin_id[i]]

# ── Row type ──────────────────────────────────────────────────────────────────

"""
    SpatialShape{G<:AbstractGeometry}

Single-shape row accessor produced by indexing into a `SpatialShapes` collection.

Carries the geometry, its `instance_id`, optional acquisition `origin`, and the
coordinate system name.

# See also
[`SpatialShapes`](@ref), [`geometry`](@ref)
"""
struct SpatialShape{G<:AbstractGeometry}
    geometry     :: G
    instance_id  :: Int32
    origin       :: Union{Nothing,String}
    coord_system :: String
end

Base.getindex(shp::SpatialShapes{G}, i::Int) where G =
    SpatialShape{G}(shp.geometries[i], shp.instance_id[i], source(shp, i), shp.coord_system)

Base.iterate(shp::SpatialShapes, i=1) = i > length(shp) ? nothing : (shp[i], i+1)
Base.eltype(::Type{SpatialShapes{G}}) where G = SpatialShape{G}

function Base.filter(pred, shp::SpatialShapes{G}) where G
    keep = [i for i in eachindex(shp.geometries) if pred(shp[i])]
    SpatialShapes(shp.geometries[keep];
                  instance_id  = shp.instance_id[keep],
                  origin_id = _subset_origin_ids(shp.origin_id, keep),
                  origin_codebook = copy(shp.origin_codebook),
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
    new_rings = map(rings) do ring
        [let v = apply(t, SVector{2,Float64}(pt[1], pt[2])); Point{2,Float32}(v[1], v[2]); end
         for pt in ring]
    end
    ext = new_rings[1]
    holes = length(new_rings) > 1 ? new_rings[2:end] : Vector{Vector{Point{2,Float32}}}()
    Polygon(ext, holes)
end

function _transform_geom(t::AbstractTransformation, multi::MultiPolygon)
    MultiPolygon([_transform_geom(t, polygon) for polygon in GeoInterface.getgeom(multi)])
end

# ── apply / apply! on SpatialPoints ──────────────────────────────────────────

function apply(t::AbstractTransformation, pts::SpatialPoints{T}) where T
    new_coords = map(pts.coords) do p
        v = apply(t, p)
        Point{2,T}(v[1], v[2])
    end
    SpatialPoints{T}(
        new_coords, copy(pts.feature_id), copy(pts.feature_codebook),
        copy(pts.instance_id), pts.feature_columns,
        isnothing(pts.origin_id) ? nothing : copy(pts.origin_id),
        copy(pts.origin_codebook), t.dst, nothing,
    )
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
    owner = _owning_dataset(pts)
    owner === nothing || touch!(owner, _dataset_ref(pts)[2])
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
    SpatialShapes(
        new_geoms;
        instance_id=copy(shp.instance_id),
        origin_id=isnothing(shp.origin_id) ? nothing : copy(shp.origin_id),
        origin_codebook=copy(shp.origin_codebook),
        coord_system=t.dst,
    )
end

function apply!(t::AbstractTransformation, shp::SpatialShapes{G}) where G
    owner = _owning_dataset(shp)
    owner === nothing || touch!(owner, _dataset_ref(shp)[2])
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
_dataset_ref(::Any)              = nothing   # immutable labels and detached extension types

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

function _clear_backref!(el::Union{SpatialPoints, SpatialShapes})
    el._attachment = nothing
    el
end
_clear_backref!(el) = el

Base.setindex!(ds::SpatialDataset, el::Union{SpatialPoints,SpatialShapes}, name::String) =
    _attach_element!(ds, el, name)

function Base.copy(pts::SpatialPoints{T}) where T
    SpatialPoints{T}(
        copy(pts.coords), copy(pts.feature_id), copy(pts.feature_codebook),
        copy(pts.instance_id), pts.feature_columns,
        isnothing(pts.origin_id) ? nothing : copy(pts.origin_id),
        copy(pts.origin_codebook), pts.coord_system, nothing,
    )
end

Base.copy(shp::SpatialShapes{G}) where G =
    SpatialShapes{G}(
        copy(shp.geometries), copy(shp.instance_id),
        isnothing(shp.origin_id) ? nothing : copy(shp.origin_id),
        copy(shp.origin_codebook), shp.coord_system, nothing,
    )

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
    SpatialPoints{T}(
        pts.coords[idx], pts.feature_id[idx], copy(pts.feature_codebook),
        pts.instance_id[idx], _subset_feature_columns(pts.feature_columns, idx),
        _subset_origin_ids(pts.origin_id, idx), copy(pts.origin_codebook),
        pts.coord_system, nothing,
    )
end

function Base.getindex(pts::SpatialPoints{T}, mask::AbstractVector{Bool}) where T
    SpatialPoints{T}(
        pts.coords[mask], pts.feature_id[mask], copy(pts.feature_codebook),
        pts.instance_id[mask], _subset_feature_columns(pts.feature_columns, mask),
        _subset_origin_ids(pts.origin_id, mask), copy(pts.origin_codebook),
        pts.coord_system, nothing,
    )
end

function Base.getindex(pts::SpatialPoints{T}, gene::String) where T
    idx = findfirst(==(gene), pts.feature_codebook)
    mask = idx === nothing ? falses(length(pts.coords)) : pts.feature_id .== Int32(idx)
    pts[mask]
end

function Base.getindex(pts::SpatialPoints, genes::AbstractVector{<:AbstractString})
    selected = Set(String.(genes))
    selected_indices = Set(
        Int32(index) for (index, gene) in pairs(pts.feature_codebook) if gene in selected
    )
    pts[BitVector(id in selected_indices for id in pts.feature_id)]
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
    count_per_instance(pts; feature=nothing) → Dict{Int32, Int}

Return a dictionary mapping each non-zero instance ID to its observation count.
When `feature` is supplied, count only observations with that feature label.

Unassigned points (`instance_id == 0`) are excluded. Useful for computing
transcript counts per cell or density metrics.

# See also
[`instance_id`](@ref), [`top_features`](@ref)
"""
function count_per_instance(pts::SpatialPoints; feature::Union{Nothing,AbstractString}=nothing)
    feature_index = if feature === nothing
        nothing
    else
        index = findfirst(==(feature), pts.feature_codebook)
        index === nothing && throw(ArgumentError(
            "feature $(repr(feature)) not found; available: $(pts.feature_codebook)",
        ))
        Int32(index)
    end
    counts = Dict{Int32, Int}()
    for index in eachindex(pts.instance_id)
        feature_index === nothing || pts.feature_id[index] == feature_index || continue
        id = pts.instance_id[index]
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
