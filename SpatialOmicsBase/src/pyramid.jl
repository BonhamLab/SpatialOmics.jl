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
    # Vector{Any} is intentional — pyramid levels are heterogeneous DiskArray
    # concrete types (ZarrV3Array, SubArray views, etc.) that share no common
    # supertype. DimensionalData and PyramidScheme.jl were evaluated and found
    # unsuitable: neither provides the Makie.Resampler callable interface with
    # coordinate normalisation needed here.
    levels::Vector{Any}
    # Global spatial extents (in physical / global-coordinate-system units).
    # Nothing → use pixel coordinates 1..nx / 1..ny (identity transform).
    x_range::Union{NTuple{2,Float64}, Nothing}
    y_range::Union{NTuple{2,Float64}, Nothing}
    # Coordinate normalisation flags (set when the image has an affine transform
    # to the global coordinate system that involves an axis swap or reflection).
    # perm_yx  : local y-axis maps to global x (and local x → global y).
    # flip_dim2: flip the post-permutation y dimension (global-y direction).
    # flip_dim3: flip the post-permutation x dimension (global-x direction).
    perm_yx::Bool
    flip_dim2::Bool
    flip_dim3::Bool
end

"""
    ImagePyramidSampler(img::SpatialImage, channel::Int=1) -> ImagePyramidSampler

Slice `channel` from each pyramid level of `img` (OME-Zarr `(c, y, x)` layout).
The slice is a lazy `view` — no pixel data is loaded until the sampler is called.

If `img.metadata` contains `"x_range"`, `"y_range"`, `"perm_yx"`, `"flip_dim2"`,
and/or `"flip_dim3"` (set by the SpatialData reader when an affine coordinate
transform is present), the sampler will render in the global coordinate system
so that overlaid images align correctly.
"""
function ImagePyramidSampler(img::SpatialImage{T}, channel::Int=1) where T
    # Channel slicing uses plain view() rather than DimensionalData — Zarr.jl
    # and PyramidScheme.jl provide no multiscale NGFF abstraction at this level.
    _slice(arr) = ndims(arr) == 3 ? view(arr, channel, :, :) :   # (ny, nx) view
                  ndims(arr) == 2 ? arr :
                  error("ImagePyramidSampler: expected 2-D or 3-D array, got $(ndims(arr))-D")
    levels = Any[_slice(img.data)]
    for lvl in img.pyramid
        push!(levels, _slice(lvl))
    end
    xr  = get(img.metadata, "x_range",   nothing)
    yr  = get(img.metadata, "y_range",   nothing)
    pyx = get(img.metadata, "perm_yx",   false)
    fd2 = get(img.metadata, "flip_dim2", false)
    fd3 = get(img.metadata, "flip_dim3", false)
    xrt = xr === nothing ? nothing : (Float64(xr[1]), Float64(xr[2]))
    yrt = yr === nothing ? nothing : (Float64(yr[1]), Float64(yr[2]))
    return ImagePyramidSampler{T}(levels, xrt, yrt, pyx, fd2, fd3)
end

# ── AbstractMatrix interface ───────────────────────────────────────────────────
# Report (nx_norm, ny_norm) so Makie creates correct (width × height) axis ranges.
# Internally levels are (ny_raw, nx_raw); if perm_yx is set the normalised
# dimensions are swapped relative to the raw storage.
function Base.size(s::ImagePyramidSampler)
    ny, nx = size(s.levels[1])   # raw (ny_raw, nx_raw)
    # When perm_yx: normalised nx = ny_raw, normalised ny = nx_raw.
    return s.perm_yx ? (ny, nx) : (nx, ny)
end

# Scalar getindex: map normalised (xi, yi) to raw (row, col) accounting for
# any axis swap or flip, then read from the finest level.
function Base.getindex(s::ImagePyramidSampler, xi::Int, yi::Int)
    finest_nx, finest_ny = size(s)
    raw_row, raw_col = if s.perm_yx
        # xi (x-norm) → raw row;  yi (y-norm) → raw col
        r = s.flip_dim3 ? (finest_nx + 1 - xi) : xi
        c = s.flip_dim2 ? (finest_ny + 1 - yi) : yi
        r, c
    else
        r = s.flip_dim2 ? (finest_ny + 1 - yi) : yi
        c = s.flip_dim3 ? (finest_nx + 1 - xi) : xi
        r, c
    end
    return s.levels[1][raw_row, raw_col]
end

# ── Callable interface for Makie.Resampler ────────────────────────────────────

"""
    (s::ImagePyramidSampler)(x::LinRange, y::LinRange) -> Matrix

Called by `Makie.Resampler` on every zoom/pan event.

`x` and `y` are **normalised pixel index coordinates** in `[1..nx_norm]` and
`[1..ny_norm]` respectively (i.e. in the same coordinate system that `size(s)`
reports).  `Makie.Resampler` always converts data / global coordinates to this
index space before calling the sampler via `resample_image`.  Code that calls
the sampler directly (e.g. `_pyramid_image!`) must do the same conversion.

Step size encodes zoom level: large step → coarser level; small step → finer.

Returns a `(length(x), length(y))` matrix — Makie places `result[i,j]` at
data coordinate `(x[i], y[j])`.
"""
function (s::ImagePyramidSampler)(x::LinRange, y::LinRange)
    finest_nx, finest_ny = size(s)   # normalised (nx_norm, ny_norm)

    # x and y are already normalised pixel indices [1..finest_N].
    # Makie.Resampler (resample_image) converts data/global coords to index
    # space before calling us; _pyramid_image! does the same for the image! path.
    xstep, ystep = step(x), step(y)

    # ── Select the best pyramid level ─────────────────────────────────────────
    best_idx  = 1
    best_dist = Inf
    for (i, lvl) in enumerate(s.levels)
        ny_k, nx_k = size(lvl)
        nx_norm_k = s.perm_yx ? ny_k : nx_k   # normalised size of this level
        ny_norm_k = s.perm_yx ? nx_k : ny_k
        dist = hypot(finest_nx / nx_norm_k - xstep, finest_ny / ny_norm_k - ystep)
        if dist < best_dist
            best_dist = dist
            best_idx  = i
        end
    end

    level      = s.levels[best_idx]
    ny_k, nx_k = size(level)
    nx_norm_k  = s.perm_yx ? ny_k : nx_k
    ny_norm_k  = s.perm_yx ? nx_k : ny_k

    # ── Normalised pixel → raw level indices ──────────────────────────────────
    # ci_norm: column indices [1..nx_norm_k] — represent the x dimension
    # ri_norm: row    indices [1..ny_norm_k] — represent the y dimension
    ci_norm = _scale_range(x, finest_nx, nx_norm_k)
    ri_norm = _scale_range(y, finest_ny, ny_norm_k)

    # Map (ci_norm, ri_norm) to raw (ri_raw, ci_raw) in the (ny_k, nx_k) level.
    ri_raw, ci_raw = if s.perm_yx
        s.flip_dim3 ? (ny_k + 1 .- ci_norm) : ci_norm,
        s.flip_dim2 ? (nx_k + 1 .- ri_norm) : ri_norm
    else
        s.flip_dim2 ? (ny_k + 1 .- ri_norm) : ri_norm,
        s.flip_dim3 ? (nx_k + 1 .- ci_norm) : ci_norm
    end

    # ── Read contiguous block from DiskArray (ascending UnitRanges) ───────────
    ri_range = minimum(ri_raw):maximum(ri_raw)
    ci_range = minimum(ci_raw):maximum(ci_raw)
    block    = collect(Float32.(level[ri_range, ci_range]))

    # ── Local indexing + orient for Makie ─────────────────────────────────────
    ri_local = ri_raw .- (first(ri_range) - 1)
    ci_local = ci_raw .- (first(ci_range) - 1)

    return s.perm_yx ? block[ri_local, ci_local] :
                       permutedims(block[ri_local, ci_local])
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
    # Display normalised pixel dimensions (accounting for perm_yx).
    szs = if s.perm_yx
        ["($(size(lvl,1))×$(size(lvl,2)))" for lvl in s.levels]
    else
        ["($(size(lvl,2))×$(size(lvl,1)))" for lvl in s.levels]
    end
    print(io, "ImagePyramidSampler{$T}($nx×$ny px, $(length(s.levels)) levels: $(join(szs, ", "))")
    s.perm_yx && print(io, ", perm_yx")
    print(io, ")")
end
