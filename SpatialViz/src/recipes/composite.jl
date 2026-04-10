# SpatialViz/src/recipes/composite.jl
#
# composite() — merge all channels of a SpatialImage into a single RGB matrix.
# spatial_panel() — multi-layer convenience figure builder.
#
# Two rendering paths:
#
#   Pyramid path (image has OME-Zarr pyramid levels)
#   ─────────────────────────────────────────────────
#   CompositePyramidSampler wraps all pyramid levels.  On each zoom/pan Makie
#   picks the right level and loads only the visible tile, then blends channels
#   to Matrix{RGBf}.  Passed to Makie.Resampler → heatmap! — same mechanism as
#   ImagePyramidSampler, which already handles Colorant element types.
#
#   Lazy path (no pyramid)
#   ─────────────────────
#   RGB mode:  ImageCore.colorview — zero-copy lazy channel stack.
#   Fluor mode: per-channel subsampled quantile pass, then MappedArrays.mappedarray
#               for lazy per-pixel blend.  Collected at render time by Makie.

using Statistics: quantile

# ---------------------------------------------------------------------------
# Shared constants / helpers
# ---------------------------------------------------------------------------

# Default colour cycle (red/green/blue first, then broader palette).
const _DEFAULT_COLORS = RGBf[
    RGBf(1.00f0, 0.00f0, 0.00f0),   # red
    RGBf(0.00f0, 0.90f0, 0.00f0),   # green
    RGBf(0.00f0, 0.00f0, 1.00f0),   # blue
    RGBf(0.00f0, 0.90f0, 0.90f0),   # cyan
    RGBf(0.90f0, 0.00f0, 0.90f0),   # magenta
    RGBf(0.90f0, 0.90f0, 0.00f0),   # yellow
]

# Resolve the `channel` + `color` kwargs to parallel (indices, colors) vectors.
#
#   channel=nothing, color=nothing  → first min(n,3) channels, red/green/blue
#   channel=nothing, color=[…]      → all channels, explicit colors (backward compat)
#   channel=[1,3]                   → those indices, default color cycle
#   channel=["DAPI","18S"]          → names resolved to indices, default colors
#   channel=[1=>:red, 3=>:cyan]     → pairs carry their own colors
#   channel=[…], color=[…]          → parallel vectors (color overrides pair colors)
#
# NOTE: inside functions that accept a `channel` kwarg the name `channel` shadows
# the exported `channels()` function — callers must use SpatialOmicsBase.channels(img).
function _resolve_channel_spec(ch_labels::Vector{String},
                                channel_spec,
                                color_arg)::Tuple{Vector{Int}, Vector{RGBf}}
    n = length(ch_labels)

    # ── no spec at all: first 3 channels, red / green / blue ──────────────────
    if channel_spec === nothing && color_arg === nothing
        n_sel = min(n, 3)
        return collect(1:n_sel), _DEFAULT_COLORS[1:n_sel]
    end

    # ── color only (backward compat): all channels, explicit colors ────────────
    if channel_spec === nothing
        len = length(color_arg)
        len == n || error("composite: $len colors for $n channels; " *
                          "pass `channel=` to select a subset")
        return collect(1:n), RGBf[RGBf(Makie.to_color(c)) for c in color_arg]
    end

    # ── channel_spec given: resolve each entry to (index, color) ──────────────
    indices    = Int[]
    pair_colors = RGBf[]
    for (j, spec) in enumerate(channel_spec)
        key, pair_color = spec isa Pair ? (spec.first, spec.second) : (spec, nothing)
        idx = if key isa Integer
            Int(key)
        elseif key isa AbstractString
            i = findfirst(==(key), ch_labels)
            i === nothing && error("composite: channel \"$key\" not found in $ch_labels")
            i
        else
            error("composite: channel key must be Int or String, got $(typeof(key))")
        end
        (1 ≤ idx ≤ n) || error("composite: channel index $idx out of range 1..$n")
        push!(indices, idx)
        push!(pair_colors, pair_color !== nothing ?
              RGBf(Makie.to_color(pair_color)) :
              _DEFAULT_COLORS[mod1(j, length(_DEFAULT_COLORS))])
    end

    # explicit color_arg overrides pair / default colors
    if color_arg !== nothing
        length(color_arg) == length(indices) ||
            error("composite: $(length(color_arg)) colors for $(length(indices)) selected channels")
        return indices, RGBf[RGBf(Makie.to_color(c)) for c in color_arg]
    end

    return indices, pair_colors
end

# Estimate the clip-th quantile from a strided subsample of a 2-D channel.
# Uses array slicing so DiskArrays loads whole chunks — fast even for Zarr-backed data.
function _subsample_quantile(ch::AbstractMatrix, clip::Float64)
    ny, nx = size(ch)
    stride = max(1, round(Int, sqrt(ny * nx / 10_000)))   # ~10 k sample points
    sample_slice = ch[1:stride:end, 1:stride:end]
    samples = vec(Float32.(collect(sample_slice)))
    isempty(samples) && return eps(Float32)
    return max(Float32(quantile(samples, clip)), eps(Float32))
end

# Mirror of _scale_range from SpatialOmicsBase/pyramid.jl — maps a LinRange in
# finest-level coordinates to integer indices in a coarser level.
function _composite_scale_range(r::LinRange, finest_N::Int, level_N::Int)::Vector{Int}
    finest_N == level_N &&
        return clamp.(round.(Int, collect(r)), 1, level_N)
    scale = (level_N - 1) / (finest_N - 1)
    return clamp.(round.(Int, (collect(r) .- 1) .* scale .+ 1), 1, level_N)
end

# ---------------------------------------------------------------------------
# CompositePyramidSampler — pyramid-aware callable for Makie.Resampler
# ---------------------------------------------------------------------------

"""
    CompositePyramidSampler{T}

Zoom-responsive composite sampler for multi-channel OME-Zarr images.

Works identically to `ImagePyramidSampler` but blends all `n` channels into
`Matrix{RGBf}` per tile, using either RGB direct mapping or fluorescence
additive blending.  Intended to be wrapped in `Makie.Resampler` and displayed
via `heatmap!`; Makie's `HeatmapShader` path handles `Colorant` element types
by skipping colormap application and rendering the RGB values directly.

Clip values for per-channel normalization are estimated from the coarsest
pyramid level at construction time — this is fast (small array) and a good
approximation for the full-resolution data.
"""
struct CompositePyramidSampler <: AbstractMatrix{RGBf}
    # All-channel levels: levels[i] is a (c, ny, nx) array, finest first.
    levels::Vector{Any}
    channel_indices::Vector{Int}   # which channels to blend (1-based into the c dim)
    colors::Vector{RGBf}
    clip_values::Vector{Float32}   # per-channel hi, estimated from coarsest level
    mode::Symbol                   # :rgb or :fluor
    # Coordinate transform fields — same semantics as ImagePyramidSampler.
    # Nothing → pixel coords 1..nx / 1..ny; set by SpatialData NGFF affine reader.
    x_range::Union{NTuple{2,Float64}, Nothing}
    y_range::Union{NTuple{2,Float64}, Nothing}
    perm_yx::Bool
    flip_dim2::Bool
    flip_dim3::Bool
end

function CompositePyramidSampler(img::SpatialImage;
                                  channel = nothing,
                                  color   = nothing,
                                  clip::Real = 0.999)
    # Use qualified name — `channel` kwarg shadows the exported channels() function.
    ch_labels = SpatialOmicsBase.channels(img)
    isempty(ch_labels) && error("CompositePyramidSampler: image has no channel axis")
    n = length(ch_labels)

    is_rgb = channel === nothing && color === nothing &&
             map(lowercase, ch_labels) == ["r", "g", "b"]

    indices, resolved = if is_rgb
        collect(1:n), RGBf[]   # color unused in RGB mode
    else
        _resolve_channel_spec(ch_labels, channel, color)
    end

    levels = Any[img.data]
    for lvl in img.pyramid
        push!(levels, lvl)
    end

    coarsest = last(levels)
    clip_values = if is_rgb
        hi = Float32(max(maximum(coarsest), eps(Float32)))
        fill(hi, n)
    else
        Float32[_subsample_quantile(view(coarsest, i, :, :), Float64(clip))
                for i in indices]
    end

    xr  = get(img.metadata, "x_range",   nothing)
    yr  = get(img.metadata, "y_range",   nothing)
    pyx = get(img.metadata, "perm_yx",   false)
    fd2 = get(img.metadata, "flip_dim2", false)
    fd3 = get(img.metadata, "flip_dim3", false)
    xrt = xr === nothing ? nothing : (Float64(xr[1]), Float64(xr[2]))
    yrt = yr === nothing ? nothing : (Float64(yr[1]), Float64(yr[2]))

    return CompositePyramidSampler(levels, indices, resolved, clip_values,
                                   is_rgb ? :rgb : :fluor,
                                   xrt, yrt, pyx, fd2, fd3)
end

# ── AbstractMatrix interface ──────────────────────────────────────────────────
# Report (nx_norm, ny_norm) — Makie heatmap convention: dim 1 → x.
# perm_yx swaps the normalised dimensions relative to raw (c, ny, nx) storage.
function Base.size(s::CompositePyramidSampler)
    _, ny, nx = size(s.levels[1])   # raw (ny_raw, nx_raw)
    return s.perm_yx ? (ny, nx) : (nx, ny)
end

# Scalar getindex — blends a single normalised-pixel from the finest level.
function Base.getindex(s::CompositePyramidSampler, xi::Int, yi::Int)
    level = s.levels[1]
    finest_nx, finest_ny = size(s)
    raw_row, raw_col = if s.perm_yx
        r = s.flip_dim3 ? (finest_nx + 1 - xi) : xi
        c = s.flip_dim2 ? (finest_ny + 1 - yi) : yi
        r, c
    else
        r = s.flip_dim2 ? (finest_ny + 1 - yi) : yi
        c = s.flip_dim3 ? (finest_nx + 1 - xi) : xi
        r, c
    end
    r = 0f0; g = 0f0; b = 0f0
    for (j, ch) in enumerate(s.channel_indices)
        v = clamp(Float32(level[ch, raw_row, raw_col]) / s.clip_values[j], 0f0, 1f0)
        col = s.colors[j]
        r += v * col.r; g += v * col.g; b += v * col.b
    end
    return RGBf(clamp(r, 0f0, 1f0), clamp(g, 0f0, 1f0), clamp(b, 0f0, 1f0))
end

# ── Callable interface for Makie.Resampler ────────────────────────────────────

"""
    (s::CompositePyramidSampler)(x::LinRange, y::LinRange) -> Matrix{RGBf}

Called by `Makie.Resampler` on every zoom/pan event. Selects the best pyramid
level, loads only the visible tile (all channels), and blends to `Matrix{RGBf}`
in Makie's `(length(x), length(y))` dimension order.
"""
function (s::CompositePyramidSampler)(x::LinRange, y::LinRange)
    # x and y are normalised pixel index coordinates [1..nx_norm] / [1..ny_norm].
    # See ImagePyramidSampler for details on the coordinate contract.
    finest_nx, finest_ny = size(s)

    xstep, ystep = step(x), step(y)

    # ── Select the best pyramid level ─────────────────────────────────────────
    best_idx  = 1
    best_dist = Inf
    for (k, lvl) in enumerate(s.levels)
        _, ny_k, nx_k = size(lvl)
        nx_norm_k = s.perm_yx ? ny_k : nx_k
        ny_norm_k = s.perm_yx ? nx_k : ny_k
        dist = hypot(finest_nx / nx_norm_k - xstep, finest_ny / ny_norm_k - ystep)
        if dist < best_dist
            best_dist = dist
            best_idx  = k
        end
    end

    level      = s.levels[best_idx]
    _, ny_k, nx_k = size(level)
    nx_norm_k  = s.perm_yx ? ny_k : nx_k
    ny_norm_k  = s.perm_yx ? nx_k : ny_k

    # ── Normalised pixel → raw level indices ──────────────────────────────────
    ci_norm = _composite_scale_range(x, finest_nx, nx_norm_k)
    ri_norm = _composite_scale_range(y, finest_ny, ny_norm_k)

    ri_raw, ci_raw = if s.perm_yx
        s.flip_dim3 ? (ny_k + 1 .- ci_norm) : ci_norm,
        s.flip_dim2 ? (nx_k + 1 .- ri_norm) : ri_norm
    else
        s.flip_dim2 ? (ny_k + 1 .- ri_norm) : ri_norm,
        s.flip_dim3 ? (nx_k + 1 .- ci_norm) : ci_norm
    end

    # ── Read contiguous block from DiskArray ───────────────────────────────────
    ri_range = minimum(ri_raw):maximum(ri_raw)
    ci_range = minimum(ci_raw):maximum(ci_raw)
    block    = collect(Float32.(level[s.channel_indices, ri_range, ci_range]))

    ri_local = ri_raw .- (first(ri_range) - 1)
    ci_local = ci_raw .- (first(ci_range) - 1)
    tile     = block[:, ri_local, ci_local]

    blend_tile = s.perm_yx ? permutedims(tile, (1, 3, 2)) : tile
    return _blend_tile(blend_tile, s.mode, s.colors, s.clip_values)
end

function Base.show(io::IO, s::CompositePyramidSampler)
    nx, ny = size(s)
    szs = if s.perm_yx
        ["($(size(lvl,2))×$(size(lvl,3)))" for lvl in s.levels]
    else
        ["($(size(lvl,3))×$(size(lvl,2)))" for lvl in s.levels]
    end
    print(io, "CompositePyramidSampler($nx×$ny px, $(length(s.levels)) levels: $(join(szs, ", ")), $(s.mode))")
    s.perm_yx && print(io, ", perm_yx")
    print(io, ")")
end

# Blend a (c, ny, nx) float tile → (nx, ny) Matrix{RGBf}.
function _blend_tile(tile::Array{Float32,3}, mode::Symbol,
                     colors::Vector{RGBf}, clip_values::Vector{Float32})
    n_ch, ny, nx = size(tile)
    out = Matrix{RGBf}(undef, nx, ny)
    if mode == :rgb
        hi = clip_values[1]
        @inbounds for y in 1:ny, x in 1:nx
            out[x, y] = RGBf(clamp(tile[1,y,x]/hi, 0f0, 1f0),
                              clamp(tile[2,y,x]/hi, 0f0, 1f0),
                              clamp(tile[3,y,x]/hi, 0f0, 1f0))
        end
    else
        cr = Float32[c.r for c in colors]
        cg = Float32[c.g for c in colors]
        cb = Float32[c.b for c in colors]
        @inbounds for y in 1:ny, x in 1:nx
            r = 0f0; g = 0f0; b = 0f0
            for i in 1:n_ch
                v = clamp(tile[i,y,x] / clip_values[i], 0f0, 1f0)
                r += v * cr[i]; g += v * cg[i]; b += v * cb[i]
            end
            out[x, y] = RGBf(clamp(r, 0f0, 1f0), clamp(g, 0f0, 1f0), clamp(b, 0f0, 1f0))
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# composite — lazy RGB matrix for non-pyramid (or explicit materialisation)
# ---------------------------------------------------------------------------

"""
    composite(img::SpatialImage; channel=nothing, color=nothing, clip=0.999) -> AbstractMatrix{RGBf}

Merge channels of `img` into a lazy `AbstractMatrix{RGBf}`.

When `img` carries OME-Zarr pyramid levels, returns a `CompositePyramidSampler`
that selects the right resolution level and loads only the visible tile on each
zoom/pan event when wrapped in `Makie.Resampler`.  The `image!` overloads do
this automatically.

For images without a pyramid (or to force materialisation), returns a lazy
`MappedArray` / `colorview` backed by the original data.  Pass the result to
`collect` to materialise, or directly to `Makie.image` / `Makie.image!`.

**Channel / colour selection** (`channel` and `color` kwargs):

- Both `nothing` (default): first 3 channels, rendered red/green/blue.
- `channel=[1,3]`: select channels by index; colours assigned from default cycle.
- `channel=["DAPI","18S"]`: select by name (must match `channels(img)` labels).
- `channel=[1=>:cyan, 3=>:magenta]`: index→colour pairs.
- `color=[:red,:blue]`: explicit colour per selected channel (overrides pair colours).

Two rendering modes selected automatically:

- **RGB mode** — channel labels exactly `["R","G","B"]` and no explicit kwargs:
  lazy `colorview`, zero allocation.
- **Fluorescence mode** — all other configs: subsampled quantile clip per channel,
  then lazy additive `mappedarray` blend.

```julia
# Pyramid-aware (zoom-responsive) display — handled automatically by image!
image!(ax, xen.images["morphology_focus"])

# First 3 channels: red/green/blue (default)
image!(ax, xen.images["morphology_focus"])

# Custom channel selection
image!(ax, img; channel=[1,3], color=[:cyan, :magenta])
image!(ax, img; channel=["DAPI", "18S"], color=[:blue, :green])
image!(ax, img; channel=[1=>:red, 2=>:green, 4=>:cyan])

# Explicit composite for a cropped region
roi = view(xen, cx-500, cx+500, cy-500, cy+500)
rgb = composite(roi.images["morphology_focus"])
image!(ax, rgb)
```
"""
function composite(img::SpatialImage;
                   channel = nothing,
                   color   = nothing,
                   clip::Real = 0.999)
    ch_labels = SpatialOmicsBase.channels(img)
    isempty(ch_labels) && error("composite: image has no channel axis")

    is_rgb = channel === nothing && color === nothing &&
             map(lowercase, ch_labels) == ["r", "g", "b"]

    if !isempty(img.pyramid)
        return CompositePyramidSampler(img; channel, color, clip)
    end

    if is_rgb
        return _composite_rgb(img.data)
    end
    indices, resolved = _resolve_channel_spec(ch_labels, channel, color)
    return _composite_fluor(img.data, indices, resolved, Float64(clip))
end

# SpatialElementView overload — crop first, then lazy-blend the region.
# Always uses the non-pyramid path since we've already bounded the region.
function composite(v::SpatialElementView{<:SpatialImage};
                   channel = nothing,
                   color   = nothing,
                   clip::Real = 0.999)
    e   = v.extent
    img = v.parent
    ax_keys = keys(img.axes)
    c_idx = findfirst(==(:c), ax_keys)
    y_idx = findfirst(==(:y), ax_keys)
    x_idx = findfirst(==(:x), ax_keys)
    (c_idx === nothing || y_idx === nothing || x_idx === nothing) &&
        error("composite: expected c/y/x axes, got $ax_keys")

    sz    = size(img.data)
    ny, nx = sz[y_idx], sz[x_idx]
    y0 = max(1, round(Int, e.ymin));  y1 = min(ny, round(Int, e.ymax))
    x0 = max(1, round(Int, e.xmin));  x1 = min(nx, round(Int, e.xmax))

    cropped = img.data[:, y0:y1, x0:x1]

    ch_labels = SpatialOmicsBase.channels(img)
    is_rgb = channel === nothing && color === nothing &&
             map(lowercase, ch_labels) == ["r", "g", "b"]

    if is_rgb
        return _composite_rgb(cropped)
    end
    indices, resolved = _resolve_channel_spec(ch_labels, channel, color)
    return _composite_fluor(cropped, indices, resolved, Float64(clip))
end

# ── internals (non-pyramid lazy path) ────────────────────────────────────────

# RGB path — lazy, zero-copy.  Normalises by global max then uses colorview.
function _composite_rgb(data::AbstractArray{<:Real, 3})
    hi = Float32(max(maximum(data), eps(Float32)))
    normed = mappedarray(x -> clamp(Float32(x) / hi, 0f0, 1f0), data)
    return colorview(RGB, view(normed, 1, :, :),
                         view(normed, 2, :, :),
                         view(normed, 3, :, :))
end

# Fluorescence path — subsampled stats pass, then lazy per-pixel blend.
# `indices` selects which channels (1-based into data's first dim) to blend.
function _composite_fluor(data::AbstractArray{<:Real, 3},
                           indices::AbstractVector{Int},
                           colors::Vector{RGBf},
                           clip::Float64)
    his = Float32[_subsample_quantile(view(data, i, :, :), clip) for i in indices]

    normed = [mappedarray(x -> clamp(Float32(x) / his[j], 0f0, 1f0),
                          view(data, i, :, :))
              for (j, i) in enumerate(indices)]

    cr = Float32[c.r for c in colors]
    cg = Float32[c.g for c in colors]
    cb = Float32[c.b for c in colors]

    blend = let cr = cr, cg = cg, cb = cb, n = length(indices)
        (vals::Vararg{Float32}) -> begin
            r = 0f0; g = 0f0; b = 0f0
            for i in Base.OneTo(n)
                v = vals[i]
                r += v * cr[i]
                g += v * cg[i]
                b += v * cb[i]
            end
            RGBf(clamp(r, 0f0, 1f0), clamp(g, 0f0, 1f0), clamp(b, 0f0, 1f0))
        end
    end

    return mappedarray(blend, normed...)
end

# ---------------------------------------------------------------------------
# image / image! overloads — SpatialImage → composite → Makie
# ---------------------------------------------------------------------------

# Pyramid-aware image display for a SpatialImage with pre-built levels.
#
# Creates an initial low-res overview via image!, then subscribes to
# ax.finallimits to update the displayed tile (xs, ys, data) whenever the
# user zooms or pans.  This bypasses Makie.Resampler / HeatmapShader, which
# would try to apply colormapping to Matrix{RGBf} and crash.
# Convert a global-coordinate viewport range to normalised pixel index space
# [1..n].  When g_min == g_max (identity transform, no x_range set) the input
# IS already in index space and is returned as-is.
function _global_to_idx(vmin::Float32, vmax::Float32,
                         g_min::Float32, g_max::Float32, n::Int)
    g_min == g_max && return vmin, vmax   # identity / no transform
    scale = (n - 1) / (g_max - g_min)
    return (clamp(1f0 + (vmin - g_min) * scale, 1f0, Float32(n)),
            clamp(1f0 + (vmax - g_min) * scale, 1f0, Float32(n)))
end

function _pyramid_image!(ax, sampler::CompositePyramidSampler; kwargs...)
    nx, ny = size(sampler)   # normalised (nx_norm, ny_norm)

    # Global spatial extents for axis positioning; sampler is always called with
    # normalised index coords [1..nx] × [1..ny] (same contract as Makie.Resampler).
    x1_full = sampler.x_range !== nothing ? Float32(sampler.x_range[1]) : 1f0
    x2_full = sampler.x_range !== nothing ? Float32(sampler.x_range[2]) : Float32(nx)
    y1_full = sampler.y_range !== nothing ? Float32(sampler.y_range[1]) : 1f0
    y2_full = sampler.y_range !== nothing ? Float32(sampler.y_range[2]) : Float32(ny)

    # Background: low-res overview in index space, displayed at global extents.
    ov_res  = min(nx, ny, 128)
    ov_tile = sampler(LinRange(1f0, Float32(nx), ov_res),
                      LinRange(1f0, Float32(ny), ov_res))
    bg = Makie.image!(ax, (x1_full, x2_full), (y1_full, y2_full), ov_tile; kwargs...)
    translate!(bg, 0, 0, -1)

    init_res  = min(nx, ny, 512)
    init_tile = sampler(LinRange(1f0, Float32(nx), init_res),
                        LinRange(1f0, Float32(ny), init_res))
    detail = Makie.image!(ax, (x1_full, x2_full), (y1_full, y2_full), init_tile; kwargs...)

    events      = ax.scene.events
    slow_limits = Observable(ax.finallimits[])

    onany(ax.finallimits, events.mousebutton, events.keyboardbutton) do lims, _, _
        isempty(events.mousebuttonstate) || return
        last = slow_limits[]
        minimum(lims) ≈ minimum(last) && widths(lims) ≈ widths(last) && return
        slow_limits[] = lims
    end

    on(slow_limits) do lims
        xmin_g, ymin_g = Float32.(minimum(lims))
        xmax_g, ymax_g = Float32.(maximum(lims))

        # Clamp viewport to image global extent.
        xmin_g = clamp(xmin_g, x1_full, x2_full);  xmax_g = clamp(xmax_g, x1_full, x2_full)
        ymin_g = clamp(ymin_g, y1_full, y2_full);  ymax_g = clamp(ymax_g, y1_full, y2_full)
        (xmax_g <= xmin_g || ymax_g <= ymin_g) && return

        # Convert global viewport to normalised pixel index space for the sampler.
        xi1, xi2 = _global_to_idx(xmin_g, xmax_g, x1_full, x2_full, nx)
        yi1, yi2 = _global_to_idx(ymin_g, ymax_g, y1_full, y2_full, ny)
        (xi2 <= xi1 || yi2 <= yi1) && return

        res_x = clamp(round(Int, xi2 - xi1), 64, 2048)
        res_y = clamp(round(Int, yi2 - yi1), 64, 2048)
        tile  = sampler(LinRange(xi1, xi2, res_x), LinRange(yi1, yi2, res_y))

        # Update displayed extents in global space; data from the index-coord tile.
        detail[1][] = Makie.EndPoints{Float32}(xmin_g, xmax_g)
        detail[2][] = Makie.EndPoints{Float32}(ymin_g, ymax_g)
        detail[3][] = tile
    end

    return detail
end

"""
    image(img; channel=nothing, color=nothing, clip=0.999, kwargs...)
    image!(ax, img; channel=nothing, color=nothing, clip=0.999, kwargs...)

Render a `SpatialImage` (or a `SpatialElementView{SpatialImage}`) as an RGB
composite. When called on a view, the axis is constrained to the view's extent.

When pyramid levels are present (OME-Zarr data), uses `_pyramid_image!` for
zoom-responsive tile loading. For images without a pyramid, calls `composite`
and collects at render time.
"""
function Makie.image(img::SpatialImageLike;
                     channel = nothing, color = nothing, clip = 0.999,
                     axis = (;), figure = (;), kwargs...)
    parent = _img_parent(img)
    fig = Figure(; figure...)
    ax  = Axis(fig[1, 1]; yreversed=true, axis...)
    if !isempty(parent.pyramid)
        sampler = CompositePyramidSampler(parent; channel, color, clip)
        plt = _pyramid_image!(ax, sampler; kwargs...)
    else
        rgb = composite(img; channel, color, clip)
        ny, nx = size(rgb)
        plt = Makie.image!(ax, 1..nx, 1..ny, permutedims(collect(rgb)); kwargs...)
    end
    _apply_extent!(ax, _img_extent(img))
    return Makie.FigureAxisPlot(fig, ax, plt)
end

function Makie.image!(ax, img::SpatialImageLike;
                      channel=nothing, color=nothing, clip=0.999, kwargs...)
    parent = _img_parent(img)
    if !isempty(parent.pyramid)
        sampler = CompositePyramidSampler(parent; channel, color, clip)
        result = _pyramid_image!(ax, sampler; kwargs...)
    else
        rgb = composite(img; channel, color, clip)
        ny, nx = size(rgb)
        result = Makie.image!(ax, 1..nx, 1..ny, permutedims(collect(rgb)); kwargs...)
    end
    _apply_extent!(ax, _img_extent(img))
    return result
end

# ---------------------------------------------------------------------------
"""
    spatial_panel(ds::SpatialDataset;
                  image_key="morphology_focus", channel=1,
                  points_key=nothing,
                  shapes_key=nothing,
                  labels_key=nothing,
                  image_colormap=:grays,
                  points_color=:red,
                  points_markersize=2,
                  shapes_color=:transparent,
                  shapes_strokecolor=:white,
                  shapes_strokewidth=0.5,
                  view::Union{SpatialDatasetView,Nothing}=nothing,
                  kwargs...)

Create a `Figure` with a single `Axis` showing the requested layers from `ds`.

Layers are rendered in order: image → labels → shapes → points, so points
always appear on top. Pass `nothing` for any key to skip that layer.

Returns `(fig, ax)`.
"""
function spatial_panel(
    ds::SpatialDataset;
    image_key::Union{String,Nothing}  = _first_key(ds.images),
    channel::Int                      = 1,
    points_key::Union{String,Nothing} = _first_key(ds.points),
    shapes_key::Union{String,Nothing} = nothing,
    labels_key::Union{String,Nothing} = nothing,
    image_colormap                    = :grays,
    points_color                      = :crimson,
    points_markersize::Real           = 2,
    shapes_color                      = :transparent,
    shapes_strokecolor                = :white,
    shapes_strokewidth::Real          = 0.5,
    view::Union{SpatialDatasetView,Nothing} = nothing,
    figure_kwargs                     = (;),
    axis_kwargs                       = (;),
    kwargs...,
)
    fig = Figure(; figure_kwargs...)
    ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true, axis_kwargs...)

    # --- image layer ---
    if image_key !== nothing && haskey(ds.images, image_key)
        img = ds.images[image_key]
        heatmap!(ax, img; channel, colormap=image_colormap)
    end

    # --- labels layer ---
    if labels_key !== nothing && haskey(ds.labels, labels_key)
        heatmap!(ax, ds.labels[labels_key]; colormap=:tab20)
    end

    # --- shapes layer ---
    if shapes_key !== nothing && haskey(ds.shapes, shapes_key)
        shp = view !== nothing ? view.shapes[shapes_key] : ds.shapes[shapes_key]
        poly!(ax, shp;
              color=shapes_color,
              strokecolor=shapes_strokecolor,
              strokewidth=shapes_strokewidth)
    end

    # --- points layer ---
    if points_key !== nothing && haskey(ds.points, points_key)
        pts = view !== nothing ? view.points[points_key] : ds.points[points_key]
        scatter!(ax, pts;
                 color=points_color,
                 markersize=points_markersize)
    end

    # Constrain axis to the view extent when one is provided; otherwise fit all layers.
    if view !== nothing
        e = view.extent
        limits!(ax, e.xmin, e.xmax, e.ymin, e.ymax)
    else
        tightlimits!(ax)
    end

    return fig, ax
end

_first_key(d::Dict) = isempty(d) ? nothing : first(keys(d))
