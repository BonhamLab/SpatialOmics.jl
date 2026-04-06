# SpatialOmicsBase/src/types/pyramid.jl
#
# ImagePyramidSampler — a callable AbstractMatrix wrapper for pre-built image
# pyramid levels (e.g. from OME-Zarr multiscales storage).
#
# Axis convention
# ───────────────
# OME-Zarr / Julia image convention:  level array has shape (height, width)
#                                     i.e. (ny, nx) = (rows, cols).
#
# Makie heatmap convention:           data[i, j] renders at (xs[i], ys[j]),
#                                     so dim 1 → x (horizontal, cols) and
#                                     dim 2 → y (vertical, rows).
#
# Therefore:
#   • levels are stored (ny, nx) — natural image layout.
#   • Base.size reports (nx, ny)  — Makie sees the correct landscape shape.
#   • The callable receives x_range (col coords) and y_range (row coords),
#     maps them to column and row indices respectively, then returns
#     permutedims(level[ri, ci])  →  (len_x, len_y) for Makie.

"""
    ImagePyramidSampler{T}

A zoom-responsive image sampler wrapping pre-built OME-Zarr pyramid levels.

Implements `AbstractMatrix{T}` with `size = (nx, ny)` (Makie heatmap convention)
and is callable as `sampler(x::LinRange, y::LinRange)` so that `Makie.Resampler`
selects the correct pyramid level on each zoom/pan event. Only the visible tile
at the chosen resolution is loaded from disk.

Levels are stored internally as `(ny, nx)` 2-D slices (one channel from the
source OME-Zarr `(c, y, x)` volume), finest first.

# Construction
```julia
sampler = ImagePyramidSampler(img)       # channel 1 (default)
sampler = ImagePyramidSampler(img, 2)    # channel 2
```
"""
struct ImagePyramidSampler{T} <: AbstractMatrix{T}
    # Stored (ny, nx): levels[1] is finest, levels[end] is coarsest.
    # Vector{Any} to accommodate heterogeneous DiskArray concrete types.
    levels::Vector{Any}
end

"""
    ImagePyramidSampler(img::SpatialImage, channel::Int=1) -> ImagePyramidSampler

Slice `channel` from each pyramid level of `img` (OME-Zarr `(c, y, x)` layout).
The slice is a lazy `view` — no pixel data is loaded until the sampler is called.
"""
function ImagePyramidSampler(img::SpatialImage{T}, channel::Int=1) where T
    _slice(arr) = ndims(arr) == 3 ? view(arr, channel, :, :) :   # (ny, nx) view
                  ndims(arr) == 2 ? arr :
                  error("ImagePyramidSampler: expected 2-D or 3-D array, got $(ndims(arr))-D")
    levels = Any[_slice(img.data)]
    for lvl in img.pyramid
        push!(levels, _slice(lvl))
    end
    return ImagePyramidSampler{T}(levels)
end

# ── AbstractMatrix interface ───────────────────────────────────────────────────
# Report (nx, ny) so Makie creates correct (width × height) axis ranges.
# Internally levels are (ny, nx), so we swap.
function Base.size(s::ImagePyramidSampler)
    ny, nx = size(s.levels[1])
    return (nx, ny)
end

# Scalar getindex: transpose access into the (ny, nx) finest level.
Base.getindex(s::ImagePyramidSampler, i::Int, j::Int) = s.levels[1][j, i]

# ── Callable interface for Makie.Resampler ────────────────────────────────────

"""
    (s::ImagePyramidSampler)(x::LinRange, y::LinRange) -> Matrix

Called by `Makie.Resampler` on every zoom/pan event.

`x` is a range of column (x / horizontal) indices in the finest-level space
`[1..nx]`; `y` is a range of row (y / vertical) indices in `[1..ny]`.

Step size encodes zoom level: large step → coarser level; small step → finer.

Returns a `(length(x), length(y))` matrix — Makie places `result[i,j]` at
data coordinate `(x[i], y[j])`.
"""
function (s::ImagePyramidSampler)(x::LinRange, y::LinRange)
    xstep, ystep = step(x), step(y)
    # finest level size in (nx, ny) terms (same convention as Base.size)
    finest_nx, finest_ny = size(s)

    # Select the pyramid level whose pixel step best matches the requested zoom.
    # For a level stored as (ny_k, nx_k):
    #   col pixel step = finest_nx / nx_k  (in finest-level x/col coordinates)
    #   row pixel step = finest_ny / ny_k  (in finest-level y/row coordinates)
    best_idx  = 1
    best_dist = Inf
    for (i, lvl) in enumerate(s.levels)
        ny_k, nx_k = size(lvl)
        dist = hypot(finest_nx / nx_k - xstep, finest_ny / ny_k - ystep)
        if dist < best_dist
            best_dist = dist
            best_idx  = i
        end
    end

    level      = s.levels[best_idx]
    ny_k, nx_k = size(level)

    # Map x (col range in finest coords [1..nx]) → column indices in this level.
    ci = _scale_range(x, finest_nx, nx_k)

    # Map y (row range in finest coords [1..ny]) → row indices in this level.
    ri = _scale_range(y, finest_ny, ny_k)

    # level[ri, ci] is (length(y), length(x)); transpose → (length(x), length(y)).
    return permutedims(level[ri, ci])
end

# Scale a LinRange from finest-level coordinates [1..finest_N] to level
# coordinates [1..level_N], returning a clamped integer index vector.
function _scale_range(r::LinRange, finest_N::Int, level_N::Int)::Vector{Int}
    if finest_N == level_N
        return clamp.(round.(Int, collect(r)), 1, level_N)
    else
        scale = (level_N - 1) / (finest_N - 1)
        return clamp.(round.(Int, (collect(r) .- 1) .* scale .+ 1), 1, level_N)
    end
end

function Base.show(io::IO, s::ImagePyramidSampler{T}) where T
    nx, ny = size(s)
    szs = ["($(size(lvl,2))×$(size(lvl,1)))" for lvl in s.levels]
    print(io, "ImagePyramidSampler{$T}($nx×$ny px, $(length(s.levels)) levels: $(join(szs, ", ")))")
end
