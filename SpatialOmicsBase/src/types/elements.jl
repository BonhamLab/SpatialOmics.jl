# SpatialOmics/src/types/elements.jl
# Abstract base and all concrete spatial element types.

"""
    SpatialElement

Abstract supertype for all spatial data elements. Every concrete element must
carry spatial coordinate information and a metadata dictionary.
"""
abstract type SpatialElement end

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

# ---------------------------------------------------------------------------
# SpatialPoints
# ---------------------------------------------------------------------------

"""
    SpatialPoints{T<:AbstractFloat}

A set of spatially-located points (e.g. transcript detections, cell centroids).

Fields
------
- `coordinates` : Matrix{T} — N×D matrix where N = num points, D = num dimensions
- `features`    : DataFrame — per-point annotations (gene, cell_id, qv, …)
- `metadata`    : Dict{String,Any}
"""
struct SpatialPoints{T<:AbstractFloat} <: SpatialElement
    coordinates::Matrix{T}
    features::DataFrame
    metadata::Dict{String,Any}
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
