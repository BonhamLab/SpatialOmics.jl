# SpatialOmicsBase/src/types/pyramid.jl
#
# ImagePyramidSampler — a callable AbstractMatrix wrapper for pre-built image
# pyramid levels (e.g. from OME-Zarr multiscales storage).
#
# Designed to be used with Makie.Resampler for zoom-responsive display:
#
#   sampler = ImagePyramidSampler(img, channel=1)
#   heatmap(x_range, y_range, Makie.Resampler(sampler))
#
# Makie.Resampler detects the callable interface via:
#   applicable(sampler, ::LinRange, ::LinRange)  →  true
# and will call sampler(x_index_range, y_index_range) on every zoom change,
# where both ranges are in the finest level's index space (1..N).
#
# Contrast with Makie.Pyramid, which builds levels eagerly from a single matrix
# using ImageBase.restrict(). ImagePyramidSampler accepts pre-computed levels
# (typically DiskArray-backed ZarrV3Arrays) and only loads the visible tile
# at the selected resolution on each render update.

"""
    ImagePyramidSampler{T}

A zoom-responsive image sampler wrapping pre-built pyramid levels.

Implements `AbstractMatrix{T}` so that Makie treats it as image data, and is
callable as `sampler(x::LinRange, y::LinRange)` so that `Makie.Resampler` uses
the pyramid path (level selection based on step size) rather than the fallback
interpolation path.

Levels are stored finest-first. Each level is a 2-D array `(ny, nx)` — a
single channel slice from the source OME-Zarr (c, y, x) volume.

# Construction
```julia
# From a SpatialImage with pyramid levels loaded
sampler = ImagePyramidSampler(img)           # channel 1
sampler = ImagePyramidSampler(img, 2)        # channel 2
```
"""
struct ImagePyramidSampler{T} <: AbstractMatrix{T}
    # Vector of 2-D arrays (ny × nx), finest [1] → coarsest [end].
    # Element type is Any to accommodate heterogeneous DiskArray concrete types.
    levels::Vector{Any}
end

"""
    ImagePyramidSampler(img::SpatialImage, channel::Int=1) -> ImagePyramidSampler

Build a sampler from a `SpatialImage`, selecting `channel` (1-based) from
the c-axis. Images are expected in OME-Zarr (c, y, x) layout; single-channel
2-D images are also accepted.
"""
function ImagePyramidSampler(img::SpatialImage{T}, channel::Int=1) where T
    function _slice(arr)
        if ndims(arr) == 3
            return view(arr, channel, :, :)   # lazy: only loads on indexing
        elseif ndims(arr) == 2
            return arr
        else
            error("ImagePyramidSampler: expected 2-D or 3-D array, got $(ndims(arr))-D")
        end
    end
    levels = Any[_slice(img.data)]
    for lvl in img.pyramid
        push!(levels, _slice(lvl))
    end
    return ImagePyramidSampler{T}(levels)
end

# AbstractMatrix interface — report the finest level's size so that Makie
# creates correct data-coordinate ranges for the heatmap.
Base.size(s::ImagePyramidSampler) = size(s.levels[1])

# Scalar indexing into the finest level (used by Makie for fallback / limits).
Base.getindex(s::ImagePyramidSampler, i::Int, j::Int) = s.levels[1][i, j]

"""
    (s::ImagePyramidSampler)(x::LinRange, y::LinRange) -> Matrix

Called by `Makie.Resampler` on every zoom change.

`x` and `y` are index ranges in the *finest-level* coordinate space (1..N).
The step sizes reflect the current zoom: large steps = zoomed out = select a
coarser level; small steps = zoomed in = select the finest level.

Returns a materialized matrix of size `(length(x), length(y))` drawn from
the best-matching pyramid level.
"""
function (s::ImagePyramidSampler)(x::LinRange, y::LinRange)
    xstep, ystep = step(x), step(y)
    finest_sz = size(s.levels[1])

    # For each level, compute the step size in finest-level coordinates.
    # Level k with size (nr_k, nc_k) has pixel step sizes:
    #   (finest_nr / nr_k,  finest_nc / nc_k)
    # Pick the level whose step size is closest to the requested (xstep, ystep).
    best_idx = 1
    best_dist = Inf
    for (i, lvl) in enumerate(s.levels)
        sz = size(lvl)
        px_step_x = finest_sz[1] / sz[1]
        px_step_y = finest_sz[2] / sz[2]
        dist = hypot(px_step_x - xstep, px_step_y - ystep)
        if dist < best_dist
            best_dist = dist
            best_idx = i
        end
    end

    level = s.levels[best_idx]
    nr, nc = size(level)
    finest_nr, finest_nc = finest_sz

    # Map from finest-level index ranges to this level's index space.
    # finest index i maps to level index: (i-1)*(nr-1)/(finest_nr-1) + 1
    if finest_nr == nr
        ri = clamp.(round.(Int, collect(x)), 1, nr)
    else
        scale_r = (nr - 1) / (finest_nr - 1)
        ri = clamp.(round.(Int, (collect(x) .- 1) .* scale_r .+ 1), 1, nr)
    end
    if finest_nc == nc
        ci = clamp.(round.(Int, collect(y)), 1, nc)
    else
        scale_c = (nc - 1) / (finest_nc - 1)
        ci = clamp.(round.(Int, (collect(y) .- 1) .* scale_c .+ 1), 1, nc)
    end

    # Materialize just this tile from the (possibly disk-backed) level array.
    return level[ri, ci]
end

function Base.show(io::IO, s::ImagePyramidSampler{T}) where T
    szs = [size(lvl) for lvl in s.levels]
    print(io, "ImagePyramidSampler{$T}($(length(s.levels)) levels: $szs)")
end
