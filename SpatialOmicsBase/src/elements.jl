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
# SpatialPoints
# ---------------------------------------------------------------------------

"""
    SpatialPoints{T<:AbstractFloat, F}

A set of spatially-located points (e.g. transcript detections, cell centroids).

`F` must satisfy the Tables.jl `istable` trait. Use a `DataFrame` for eager
in-memory data, or a lazy handle (e.g. `Parquet2.Dataset`) to defer loading
until the data are actually accessed.

Fields
------
- `coordinates` : AbstractMatrix{T} — N×D matrix where N = num points, D = num dimensions
- `features`    : F — per-point annotations (gene, cell_id, qv, …); Tables.istable
- `metadata`    : Dict{String,Any}
"""
struct SpatialPoints{T<:AbstractFloat, F} <: SpatialElement
    coordinates::AbstractMatrix{T}
    features::F
    metadata::Dict{String,Any}

    function SpatialPoints(coords::AbstractMatrix{<:AbstractFloat}, features, meta::Dict{String,Any})
        Tables.istable(features) ||
            throw(ArgumentError("SpatialPoints: features must satisfy Tables.istable"))
        T = eltype(coords)
        new{T, typeof(features)}(coords, features, meta)
    end
end

# Convenience constructor that infers T from element type
SpatialPoints(coords, features, meta) =
    SpatialPoints(coords, features, convert(Dict{String,Any}, meta))

function Base.show(io::IO, pts::SpatialPoints{T}) where T
    n, d = size(pts.coordinates)
    cols = Tables.columnnames(pts.features)
    ncols = length(cols)
    limit = 4
    col_str = ncols <= limit ? join(cols, ", ") :
              join(cols[1:limit-1], ", ") * ", …+$(ncols-(limit-1))"
    print(io, "SpatialPoints{$T}($n × $(d)D, features: [$col_str])")
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
# SpatialShapes
# ---------------------------------------------------------------------------

"""
    SpatialShapes

A collection of geometric shapes (circles, polygons) representing cell
boundaries or tissue regions as vectors, rather than rasterised label images.

Fields
------
- `geometries` : Vector — each element is a shape description (polygon / circle)
- `features`   : DataFrame — per-shape annotations
- `metadata`   : Dict{String,Any}
"""
struct SpatialShapes <: SpatialElement
    geometries::Vector{Any}
    features::DataFrame
    metadata::Dict{String,Any}
end

function Base.show(io::IO, shp::SpatialShapes)
    print(io, "SpatialShapes(", length(shp.geometries), " shapes, features: [", _cols_str(shp.features), "])")
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
