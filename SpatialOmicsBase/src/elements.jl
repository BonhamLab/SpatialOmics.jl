# SpatialOmics/src/types/elements.jl
# Abstract base and all concrete spatial element types.

"""
    SpatialElement

Abstract supertype for all spatial data elements. Every concrete element must
carry spatial coordinate information and a metadata dictionary.
"""
abstract type SpatialElement end

# Shared display helper — truncates a DataFrame's column list to avoid long lines.
_cols_str(df::DataFrame, limit::Int=4) =
    ncol(df) <= limit ? join(names(df), ", ") :
                        join(names(df)[1:limit-1], ", ") * ", …+$(ncol(df)-(limit-1))"

# ---------------------------------------------------------------------------
# SpatialImage
# ---------------------------------------------------------------------------

"""
    SpatialImage{T<:AbstractFloat}

A spatially-resolved image (e.g. DAPI, immunofluorescence, H&E).
`data` supports lazy DiskArray / memory-mapped arrays so large images are
never fully loaded into RAM unless explicitly requested.

OME-Zarr images are stored in `(c, y, x)` axis order. Multi-resolution stores
expose additional coarser levels via `pyramid` (finest = `data`, then coarser
in `pyramid[1]`, `pyramid[2]`, …). Use `ImagePyramidSampler(img, channel)` to
build a zoom-responsive sampler for `Makie.Resampler`.

Fields
------
- `data`     : AbstractArray{T, N} — full-resolution pixel data (supports DiskArray)
- `pyramid`  : Vector{Any} — coarser pyramid levels; empty if single-resolution
- `axes`     : NamedTuple — spatial axes with units, e.g. `(c=..., y=..., x=...)`
- `metadata` : Dict{String,Any}
"""
struct SpatialImage{T} <: SpatialElement
    data::AbstractArray{T}
    pyramid::Vector{Any}   # coarser-resolution levels, finest-excluded; may be DiskArrays
    axes::NamedTuple
    metadata::Dict{String,Any}
end

# Backward-compatible 3-arg constructor (no pyramid).
SpatialImage(data::AbstractArray{T}, axes::NamedTuple, metadata::Dict{String,Any}) where T =
    SpatialImage{T}(data, Any[], axes, metadata)

"""
    SpatialImage(data::AbstractArray{T,3}, channels::AbstractVector{<:AbstractString})

Construct a `SpatialImage` from a `(c, y, x)` array and channel name vector.
Synthesises standard `(c, y, x)` axes; stores `channels` in metadata.
"""
function SpatialImage(data::AbstractArray{T,3},
                      channels::AbstractVector{<:AbstractString}) where T
    size(data, 1) == length(channels) || throw(ArgumentError(
        "SpatialImage: $(length(channels)) channel names for $(size(data, 1)) channels"))
    SpatialImage{T}(data, Any[],
                    (c = nothing, y = nothing, x = nothing),
                    Dict{String,Any}("channel_labels" => collect(String, channels)))
end

"""
    SpatialImage(img::AbstractMatrix{C}) where C <: Colorant

Construct a `SpatialImage` from a color pixel matrix (e.g. `FileIO.load("image.jpg")`).
Uses `reinterpret(reshape, ...)` to produce a zero-copy `(c, y, x)` channel array.
Channel names default to the color component field names (e.g. `["r", "g", "b"]` for RGB).
"""
function SpatialImage(img::AbstractMatrix{C}) where {C <: Colorant}
    T  = eltype(C)
    fnames = fieldnames(C)
    nc = length(fnames)
    h, w = size(img)
    data = nc == 1 ? reshape(reinterpret(T, img), 1, h, w) :
                     reinterpret(reshape, T, img)
    channels = nc == 1 ? [string(nameof(C))] : [string(f) for f in fnames]
    SpatialImage(collect(data), channels)
end

"""
    channels(img::SpatialImage) -> Vector{String}

Return the channel labels for `img`.

Labels are read from the `omero.channels[i].label` section of the OME-NGFF
metadata when present. Falls back to `["ch_1", "ch_2", …]` if the metadata is
absent or incomplete. Returns an empty vector for images with no channel axis
(i.e. `img.axes` has no `:c` key).
"""
function channels(img::SpatialImage)
    c_idx = findfirst(==(:c), keys(img.axes))
    c_idx === nothing && return String[]
    n = size(img.data, c_idx)

    # User-set labels take priority over anything in the zarr metadata
    user = get(img.metadata, "channel_labels", nothing)
    user !== nothing && length(user) == n && return convert(Vector{String}, user)

    lbls = _omero_channel_labels(img.metadata, n)
    lbls !== nothing && return lbls

    return ["ch_$i" for i in 1:n]
end

"""
    channels!(img::SpatialImage, labels::Vector{String}) -> SpatialImage

Set the channel labels for `img`. Stored in `img.metadata["channel_labels"]`
and returned by `channels(img)` with higher priority than OME-NGFF metadata.

```julia
channels!(images(xen, "he_image"), ["R", "G", "B"])
```
"""
function channels!(img::SpatialImage, labels::Vector{String})
    c_idx = findfirst(==(:c), keys(img.axes))
    c_idx === nothing && error("channels!: image has no channel axis")
    n = size(img.data, c_idx)
    length(labels) == n ||
        error("channels!: $(length(labels)) labels for $n channels")
    img.metadata["channel_labels"] = labels
    return img
end

# Walk the zarr_attrs → attributes → omero → channels path without try/catch.
# Returns a Vector{String} of length n, or nothing if the metadata is absent/incomplete.
function _omero_channel_labels(meta::Dict{String,Any}, n::Int)
    zarr_attrs = get(meta, "zarr_attrs", nothing)
    zarr_attrs === nothing && return nothing
    haskey(zarr_attrs, :attributes) || return nothing
    attrs = zarr_attrs.attributes
    # OME-NGFF stores omero under attributes.ome.omero; fall back to
    # attributes.omero for files that omit the intermediate :ome key.
    omero = if haskey(attrs, :ome) && haskey(attrs.ome, :omero)
        attrs.ome.omero
    elseif haskey(attrs, :omero)
        attrs.omero
    else
        return nothing
    end
    haskey(omero, :channels) || return nothing
    ch = omero.channels
    length(ch) == n || return nothing
    haskey(first(ch), :label) || return nothing
    return [String(c.label) for c in ch]
end

function Base.show(io::IO, img::SpatialImage{T}) where T
    sz = size(img.data)
    dims = if !isempty(img.axes)
        join(["$k=$(sz[i])" for (i, k) in enumerate(keys(img.axes))], " × ")
    else
        join(sz, "×")
    end
    nlevels = length(img.pyramid)
    pyramid_str = nlevels > 0 ? ", $nlevels pyramid level$(nlevels == 1 ? "" : "s")" : ""
    ch = channels(img)
    ch_str = isempty(ch) ? "" : ", channels: [$(join(ch, ", "))]"
    print(io, "SpatialImage{$T}($dims$pyramid_str$ch_str)")
end

# ---------------------------------------------------------------------------
# SpatialPoints / SpatialPoint
# ---------------------------------------------------------------------------

"""
    SpatialPoint{T<:AbstractFloat}

Single spatially-located point. Row-accessor returned by `getindex` on
`SpatialPoints`; not typically constructed directly.

Fields
------
- `coordinates` : NTuple{2,T} — (x, y)
- `label`       : String — per-point label (e.g. gene name); empty string if none
- `features`    : Union{Nothing,NamedTuple} — extra annotations for this row
"""
struct SpatialPoint{T<:AbstractFloat}
    coordinates :: NTuple{2,T}
    label       :: String
    features    :: Union{Nothing, NamedTuple}
end

"""
    SpatialPoints{T<:AbstractFloat} <: SpatialElement

Set of spatially-located points (e.g. transcript detections, cell centroids).

Storage is struct-of-arrays: `coordinates` is N×2 `Matrix{T}`, `labels` is
`Vector{String}` (empty if no labels), and `index` maps each unique label to its
row indices for O(1) lookup via `pts["gene"]` or `pts[["g1","g2"]]`.

`features` holds extra columns beyond x/y/label. The package never reads it
internally — it is purely for user access. Views can index into it if not
`nothing`.

# Future optimization
`features` is eagerly materialised. A lazy path could store a Tables.jl source
(e.g. `Parquet2.Dataset`) directly; the primary constructor is the natural place
to add that branch.

Fields
------
- `coordinates` : Matrix{T} — N×2
- `labels`      : Vector{String} — per-point label; may be empty
- `features`    : Union{Nothing,NamedTuple} — extra columns; not used internally
- `index`       : Dict{String,Vector{Int}} — label → row indices; built from labels
- `metadata`    : Dict{String,Any}
"""
struct SpatialPoints{T<:AbstractFloat} <: SpatialElement
    coordinates :: Matrix{T}
    labels      :: Vector{String}
    features    :: Union{Nothing, NamedTuple}
    index       :: Dict{String, Vector{Int}}
    metadata    :: Dict{String,Any}
end

function _build_label_index(labels::Vector{String})::Dict{String,Vector{Int}}
    idx = Dict{String,Vector{Int}}()
    for (i, v) in enumerate(labels)
        push!(get!(idx, v, Int[]), i)
    end
    return idx
end

# Internal constructor — always rebuilds index from labels.
function SpatialPoints{T}(coords::Matrix{T},
                           labels::Vector{String},
                           features::Union{Nothing,NamedTuple},
                           meta::Dict{String,Any}) where T<:AbstractFloat
    SpatialPoints{T}(coords, labels, features, _build_label_index(labels), meta)
end

"""
    SpatialPoints(table; x_col=:x, y_col=:y, label_col=nothing, meta=Dict{String,Any}())

Construct a `SpatialPoints` from any Tables.jl source. `x_col`/`y_col` become
the coordinate matrix; `label_col` (if given) becomes `labels` and the O(1)
index; all remaining columns become `features`.

Column name arguments accept `Symbol` or `String`.

```julia
SpatialPoints(df; x_col=:x_global_px, y_col=:y_global_px, label_col=:gene)
SpatialPoints(df)                    # defaults to :x, :y; no labels
```
"""
function SpatialPoints(table;
                       x_col     = :x,
                       y_col     = :y,
                       label_col = nothing,
                       meta      :: Dict{String,Any} = Dict{String,Any}())
    Tables.istable(table) ||
        throw(ArgumentError("SpatialPoints: table must satisfy Tables.istable"))
    x_sym = Symbol(x_col)
    y_sym = Symbol(y_col)
    xs = Tables.getcolumn(table, x_sym)
    ys = Tables.getcolumn(table, y_sym)
    T  = promote_type(eltype(xs), eltype(ys), Float32)
    coords = Matrix{T}(undef, length(xs), 2)
    coords[:, 1] .= xs
    coords[:, 2] .= ys
    lbl_sym = label_col === nothing ? nothing : Symbol(label_col)
    labels  = lbl_sym === nothing ? String[] :
              Vector{String}(string.(Tables.getcolumn(table, lbl_sym)))
    skip = lbl_sym === nothing ? (x_sym, y_sym) : (x_sym, y_sym, lbl_sym)
    feat_names = Tuple(Symbol(nm) for nm in Tables.columnnames(table) if Symbol(nm) ∉ skip)
    features   = isempty(feat_names) ? nothing :
                 NamedTuple{feat_names}(Tuple(Tables.getcolumn(table, nm) for nm in feat_names))
    SpatialPoints{T}(coords, labels, features, meta)
end

"""
    SpatialPoints(points::AbstractVector{<:GeometryBasics.Point2}; meta=Dict{String,Any}())

Construct from a vector of `GeometryBasics.Point2`. No labels, no features.
"""
function SpatialPoints(points::AbstractVector{<:GeometryBasics.Point2};
                       meta::Dict{String,Any} = Dict{String,Any}())
    T      = eltype(eltype(points))
    coords = Matrix{T}(undef, length(points), 2)
    for (i, p) in enumerate(points)
        coords[i, 1] = p[1]
        coords[i, 2] = p[2]
    end
    SpatialPoints{T}(coords, String[], nothing, meta)
end

Base.length(pts::SpatialPoints)    = size(pts.coordinates, 1)
Base.eachindex(pts::SpatialPoints) = 1:length(pts)

"""
    labels(pts::SpatialPoints) -> Vector{String}

Return the per-point label vector (e.g. gene names). Empty if no labels were
set at construction. Mirrors `labels(ds::SpatialDataset)` semantics.
"""
labels(pts::SpatialPoints) = pts.labels

# Row accessor
function Base.getindex(pts::SpatialPoints{T}, i::Int) where T
    feats = pts.features === nothing ? nothing :
            NamedTuple{keys(pts.features)}(Tuple(v[i] for v in values(pts.features)))
    SpatialPoint{T}((pts.coordinates[i,1], pts.coordinates[i,2]),
                    isempty(pts.labels) ? "" : pts.labels[i],
                    feats)
end

# Bool/BitVector mask — rebuilds index
function Base.getindex(pts::SpatialPoints{T}, mask::Union{AbstractVector{Bool},BitVector}) where T
    feats = pts.features === nothing ? nothing :
            NamedTuple{keys(pts.features)}(Tuple(v[mask] for v in values(pts.features)))
    SpatialPoints{T}(pts.coordinates[mask, :],
                     isempty(pts.labels) ? String[] : pts.labels[mask],
                     feats, pts.metadata)
end

# Integer vector — rebuilds index
function Base.getindex(pts::SpatialPoints{T}, indices::AbstractVector{<:Integer}) where T
    feats = pts.features === nothing ? nothing :
            NamedTuple{keys(pts.features)}(Tuple(v[indices] for v in values(pts.features)))
    SpatialPoints{T}(pts.coordinates[indices, :],
                     isempty(pts.labels) ? String[] : pts.labels[indices],
                     feats, pts.metadata)
end

# String label — O(1) via index
Base.getindex(pts::SpatialPoints, key::AbstractString) =
    pts[get(pts.index, key, Int[])]

# String vector — O(1) multi-label union
Base.getindex(pts::SpatialPoints, keys::AbstractVector{<:AbstractString}) =
    pts[sort!(mapreduce(k -> get(pts.index, k, Int[]), vcat, keys))]

Base.iterate(pts::SpatialPoints, i=1) = i > length(pts) ? nothing : (pts[i], i+1)

"""
    Base.filter(f, pts::SpatialPoints) -> SpatialPoints

Return subset of `pts` for which predicate `f(::SpatialPoint)` is true.

```julia
filter(pt -> pt.label == "Spry2", pts)
filter(pt -> pt.features.cell_id ∈ valid_ids, pts)
```
"""
function Base.filter(f, pts::SpatialPoints)
    mask = [f(pts[i]) for i in 1:length(pts)]
    pts[mask]
end

function Base.show(io::IO, pts::SpatialPoints{T}) where T
    n      = size(pts.coordinates, 1)
    nlbls  = length(pts.index)
    feats  = pts.features === nothing ? "" : begin
        cols  = keys(pts.features)
        ncols = length(cols)
        limit = 4
        s = ncols <= limit ? join(cols, ", ") :
            join(cols[1:limit-1], ", ") * ", …+$(ncols-(limit-1))"
        ", features: [$s]"
    end
    lbl_str = isempty(pts.labels) ? "" : ", $nlbls unique labels"
    print(io, "SpatialPoints{$T}($n × 2D$lbl_str$feats)")
end

# ---------------------------------------------------------------------------
# SpatialLabels
# ---------------------------------------------------------------------------

"""
    SpatialLabels

An integer label image where each unique integer identifies a segmentation
instance (cell, nucleus, tissue region).

Fields
------
- `data`     : AbstractArray{Int} — same spatial extent as the parent image
- `metadata` : Dict{String,Any}
"""
struct SpatialLabels <: SpatialElement
    data::AbstractArray{<:Integer}
    metadata::Dict{String,Any}
end

function Base.show(io::IO, lbl::SpatialLabels)
    print(io, "SpatialLabels(", join(size(lbl.data), "×"), " ", eltype(lbl.data), ")")
end

# ---------------------------------------------------------------------------
# SpatialShape / SpatialShapes
# ---------------------------------------------------------------------------

"""
    SpatialShape{G}

A single geometric shape with a string name and optional metadata.
`G` must be a `GeometryBasics.AbstractGeometry` subtype or `Nothing`
(used when geometry is unavailable, e.g. after an HDF5 round-trip that
stores only tabular data).

Fields
------
- `geometry` : G — the geometry (`GeometryBasics.Polygon`, `GeometryBasics.Circle`, `nothing`, …)
- `name`     : String — primary identifier (cell name, ROI label, FOV id, …)
- `metadata` : NamedTuple — per-shape annotations; schema may vary across shapes

Convenience constructors set `name` to a random 8-character string and
`metadata` to an empty `NamedTuple` when omitted:

```julia
SpatialShape(poly)                          # random name, no metadata
SpatialShape(poly, "roi_ventricle")         # explicit name, no metadata
SpatialShape(poly, "cell_42", (area=12.3,)) # full explicit form
```
"""
struct SpatialShape{G <: Union{Nothing, GeometryBasics.AbstractGeometry}}
    geometry::G
    name::String
    metadata::NamedTuple
end

const _SHAPE_NAME_CHARS = collect("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
_random_shape_name() = String(rand(_SHAPE_NAME_CHARS, 8))

SpatialShape(g::G) where {G <: Union{Nothing, GeometryBasics.AbstractGeometry}} =
    SpatialShape{G}(g, _random_shape_name(), NamedTuple())

SpatialShape(g::G, name::String) where {G <: Union{Nothing, GeometryBasics.AbstractGeometry}} =
    SpatialShape{G}(g, name, NamedTuple())

"""
    SpatialShapes{G}

A collection of geometric shapes sharing geometry type `G`.

Metadata schemas may differ across shapes (heterogeneous `NamedTuple`s are
widened to a common column set when the Tables.jl interface is used).

Implements Tables.jl, so `DataFrame(ss)` works automatically. The resulting
table always has a `name` column followed by the union of all metadata keys.

Fields
------
- `shapes`   : Vector{SpatialShape{G}}
- `metadata` : Dict{String,Any}
"""
struct SpatialShapes{G <: Union{Nothing, GeometryBasics.AbstractGeometry}} <: SpatialElement
    shapes::Vector{SpatialShape{G}}
    metadata::Dict{String,Any}
end

SpatialShapes(shapes::Vector{SpatialShape{G}}) where {G} =
    SpatialShapes{G}(shapes, Dict{String,Any}())

# Preferred name columns when building SpatialShapes from a tabular features argument.
# String-typed identifiers (name, cell) take priority over numeric IDs (cell_id).
const _SHAPE_NAME_COLS = (:name, :cell, :cell_id, :fov_id, :id)

"""
    SpatialShapes(geoms, features, metadata)

Backward-compatible constructor. Pairs each geometry with the corresponding row
of `features` (any Tables.jl source). The first column whose name appears in
`$_SHAPE_NAME_COLS` becomes `shape.name` (converted to `String`). All feature
columns are also stored verbatim in `shape.metadata` so downstream join keys
(e.g. `instance_key`) continue to work against their original types.
"""
function SpatialShapes(geoms::AbstractVector{G}, features, metadata::Dict{String,Any}) where {G <: Union{Nothing, GeometryBasics.AbstractGeometry}}
    rows = collect(Tables.rows(features))
    isempty(rows) && return SpatialShapes{G}(SpatialShape{G}[], metadata)
    col_syms = Tuple(Tables.columnnames(first(rows)))
    nk_idx   = findfirst(k -> k ∈ col_syms, _SHAPE_NAME_COLS)
    name_key = nk_idx !== nothing ? _SHAPE_NAME_COLS[nk_idx] : nothing
    shapes = [begin
        nm   = name_key !== nothing ? string(Tables.getcolumn(r, name_key)) : ""
        meta = NamedTuple{col_syms}(Tuple(Tables.getcolumn(r, c) for c in col_syms))
        SpatialShape{G}(g, nm, meta)
    end for (g, r) in zip(geoms, rows)]
    return SpatialShapes{G}(shapes, metadata)
end

# ---------------------------------------------------------------------------
# Tables.jl interface — enables DataFrame(ss)
# ---------------------------------------------------------------------------

# Union of all metadata keys across shapes.  Fast-paths the common case where
# all shapes share the same schema.
function _meta_keys(ss::SpatialShapes)
    isempty(ss.shapes) && return ()
    fk = keys(first(ss.shapes).metadata)
    all(s -> keys(s.metadata) === fk, ss.shapes) && return fk
    kset = Set{Symbol}(fk)
    for s in @view ss.shapes[2:end]
        union!(kset, keys(s.metadata))
    end
    return Tuple(sort!(collect(kset)))
end

Tables.istable(::Type{<:SpatialShapes}) = true
Tables.columnaccess(::Type{<:SpatialShapes}) = true

function Tables.columns(ss::SpatialShapes)
    name_col = [s.name for s in ss.shapes]
    mk = _meta_keys(ss)
    isempty(mk) && return (name = name_col,)
    meta_cols = NamedTuple{mk}(
        Tuple([get(s.metadata, k, missing) for s in ss.shapes] for k in mk))
    return merge((name = name_col,), meta_cols)
end

Tables.columnnames(ss::SpatialShapes) = (:name, _meta_keys(ss)...)

# `name` is a first-class column; other names delegate to per-shape metadata.
Tables.getcolumn(ss::SpatialShapes, nm::Symbol) =
    nm === :name ? [s.name for s in ss.shapes] :
                   [get(s.metadata, nm, missing) for s in ss.shapes]
Tables.getcolumn(ss::SpatialShapes, i::Int) = getfield(Tables.columns(ss), i)

# Collection interface
Base.length(ss::SpatialShapes)           = length(ss.shapes)
Base.getindex(ss::SpatialShapes, i::Int) = ss.shapes[i]
Base.iterate(ss::SpatialShapes, args...) = iterate(ss.shapes, args...)

"""
    geometry(ss::SpatialShapes, i::Int) -> G

Return the geometry of shape `i` in `ss`.
"""
geometry(ss::SpatialShapes, i::Int) = ss.shapes[i].geometry

"""
    geometries(ss::SpatialShapes) -> Vector

Return all geometries in `ss`.
"""
geometries(ss::SpatialShapes) = [s.geometry for s in ss.shapes]

function Base.show(io::IO, shp::SpatialShapes{G}) where G
    mk = _meta_keys(shp)
    meta_str = isempty(mk) ? "" :
               length(mk) <= 3 ? ", meta: [$(join(string.(mk), ", "))]" :
               ", meta: [$(join(string.(mk[1:3]), ", ")), …+$(length(mk)-3)]"
    print(io, "SpatialShapes{$G}($(length(shp.shapes)) shapes$meta_str)")
end

"""
    Base.in(pt::SpatialPoint, roi::SpatialShapes) -> Bool

Test whether `pt` lies inside any shape in `roi` using signed distance.
Enables `pt ∈ roi` and `filter(pt -> pt ∈ roi, pts)`.
"""
function Base.in(pt::SpatialPoint, roi::SpatialShapes)
    p = GeometryBasics.Point2(pt.coordinates...)
    any(s.geometry !== nothing &&
        GeometryOps.signed_distance(p, s.geometry) <= 0
        for s in roi.shapes)
end

# ---------------------------------------------------------------------------
# SpatialTable
# ---------------------------------------------------------------------------

"""
    SpatialTable

A feature-by-observation table (analogous to AnnData's X + obs + var).
Typically stores gene expression counts aligned to spatial observations.

Fields
------
- `data`     : AbstractMatrix — observations × features count matrix
- `obs`      : DataFrame — per-observation (cell) annotations
- `var`      : DataFrame — per-feature (gene) annotations
- `metadata` : Dict{String,Any}
"""
struct SpatialTable <: SpatialElement
    data::AbstractMatrix
    obs::DataFrame
    var::DataFrame
    metadata::Dict{String,Any}
end

function Base.show(io::IO, tbl::SpatialTable)
    n_obs, n_var = size(tbl.data)
    print(io, "SpatialTable($n_obs obs × $n_var vars, obs: [", _cols_str(tbl.obs), "])")
end

# ---------------------------------------------------------------------------
# SpatialTable metadata accessors
# ---------------------------------------------------------------------------

# Read a scalar string from zarr_attrs.attributes, returning `default` if absent.
function _spatialdata_attr(meta::Dict{String,Any}, key::Symbol, default)
    zarr_attrs = get(meta, "zarr_attrs", nothing)
    zarr_attrs === nothing && return default
    haskey(zarr_attrs, :attributes) || return default
    attrs = zarr_attrs.attributes
    haskey(attrs, key) || return default
    val = attrs[key]
    return val isa AbstractString ? String(val) : val
end

"""
    region(tbl::SpatialTable) -> Union{String, Nothing}

Return the name of the spatial element this table annotates — the `region`
field from the SpatialData NGFF metadata (e.g. `"cell_circles"`).
Returns `nothing` if the metadata is absent.
"""
region(tbl::SpatialTable) = _spatialdata_attr(tbl.metadata, :region, nothing)

"""
    region_key(tbl::SpatialTable) -> String

Return the `obs` column name that identifies which region each row belongs to.
Defaults to `"region"` when metadata is absent.
"""
region_key(tbl::SpatialTable) = _spatialdata_attr(tbl.metadata, :region_key, "region")

"""
    instance_key(tbl::SpatialTable) -> String

Return the `obs` column name that links each row to a specific instance in
the linked region element (e.g. `"cell_id"`).
Defaults to `"instance_id"` when metadata is absent.
"""
instance_key(tbl::SpatialTable) = _spatialdata_attr(tbl.metadata, :instance_key, "instance_id")
