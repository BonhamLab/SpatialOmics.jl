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

# Default per-channel colours for fluorescence (cycles if more channels than entries).
const _FLUOR_DEFAULTS = RGBf[
    RGBf(0.20f0, 0.40f0, 1.00f0),   # blue   — DAPI / nuclear stains
    RGBf(0.10f0, 0.90f0, 0.10f0),   # green  — FITC / GFP
    RGBf(1.00f0, 0.15f0, 0.15f0),   # red    — TRITC / mCherry
    RGBf(0.00f0, 0.90f0, 0.90f0),   # cyan
    RGBf(0.90f0, 0.00f0, 0.90f0),   # magenta
    RGBf(0.90f0, 0.90f0, 0.00f0),   # yellow
]

function _resolve_composite_colors(colors, n::Int)::Vector{RGBf}
    colors === nothing &&
        return [_FLUOR_DEFAULTS[mod1(i, length(_FLUOR_DEFAULTS))] for i in 1:n]
    length(colors) == n ||
        error("composite: $(length(colors)) colours for $n channels")
    return [RGBf(Makie.to_color(c)) for c in colors]
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
    colors::Vector{RGBf}
    clip_values::Vector{Float32}   # per-channel hi, estimated from coarsest level
    mode::Symbol                   # :rgb or :fluor
end

function CompositePyramidSampler(img::SpatialImage;
                                  colors = nothing,
                                  clip::Real = 0.999)
    ch_labels = channels(img)
    isempty(ch_labels) && error("CompositePyramidSampler: image has no channel axis")
    n = length(ch_labels)

    is_rgb   = colors === nothing && map(lowercase, ch_labels) == ["r", "g", "b"]
    resolved = _resolve_composite_colors(colors, n)

    # Build level list: base data first, then coarser pyramid levels.
    # All kept in full (c, ny, nx) layout — channel slicing happens per tile.
    levels = Any[img.data]
    for lvl in img.pyramid
        push!(levels, lvl)
    end

    # Estimate clip values from the coarsest (smallest, fastest-to-read) level.
    coarsest = last(levels)
    clip_values = if is_rgb
        # RGB: single global max applied to all channels to preserve colour balance.
        hi = Float32(max(maximum(coarsest), eps(Float32)))
        fill(hi, n)
    else
        Float32[_subsample_quantile(view(coarsest, i, :, :), Float64(clip)) for i in 1:n]
    end

    return CompositePyramidSampler(levels, resolved, clip_values, is_rgb ? :rgb : :fluor)
end

# ── AbstractMatrix interface ──────────────────────────────────────────────────
# Report (nx, ny) — Makie heatmap convention: dim 1 → x.
function Base.size(s::CompositePyramidSampler)
    _, ny, nx = size(s.levels[1])
    return (nx, ny)
end

# Scalar getindex — blends a single pixel from the finest level.
function Base.getindex(s::CompositePyramidSampler, xi::Int, yi::Int)
    level  = s.levels[1]
    n_ch   = length(s.colors)
    r = 0f0; g = 0f0; b = 0f0
    for i in 1:n_ch
        v = clamp(Float32(level[i, yi, xi]) / s.clip_values[i], 0f0, 1f0)
        c = s.colors[i]
        r += v * c.r; g += v * c.g; b += v * c.b
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
    finest_nx, finest_ny = size(s)

    # Select the pyramid level whose pixel step best matches the requested zoom.
    best_idx  = 1
    best_dist = Inf
    for (k, lvl) in enumerate(s.levels)
        _, ny_k, nx_k = size(lvl)
        dist = hypot(finest_nx / nx_k - step(x), finest_ny / ny_k - step(y))
        if dist < best_dist
            best_dist = dist
            best_idx  = k
        end
    end

    level = s.levels[best_idx]
    _, ny_k, nx_k = size(level)

    ci = _composite_scale_range(x, finest_nx, nx_k)
    ri = _composite_scale_range(y, finest_ny, ny_k)

    ri_range = minimum(ri):maximum(ri)
    ci_range = minimum(ci):maximum(ci)

    # Load all channels for this tile as a concrete (c, ny_tile, nx_tile) array.
    block = collect(Float32.(level[:, ri_range, ci_range]))

    ri_local = ri .- (first(ri_range) - 1)
    ci_local = ci .- (first(ci_range) - 1)
    tile = block[:, ri_local, ci_local]   # (c, ny_out, nx_out)

    # Blend to (nx_out, ny_out) Matrix{RGBf} — Makie's data[i,j] @ (x[i], y[j]).
    return _blend_tile(tile, s.mode, s.colors, s.clip_values)
end

function Base.show(io::IO, s::CompositePyramidSampler)
    nx, ny = size(s)
    szs = ["($(size(lvl,3))×$(size(lvl,2)))" for lvl in s.levels]
    print(io, "CompositePyramidSampler($nx×$ny px, $(length(s.levels)) levels: $(join(szs, ", ")), $(s.mode))")
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
    composite(img::SpatialImage; colors=nothing, clip=0.999) -> AbstractMatrix{RGBf}

Merge all channels of `img` into a lazy `AbstractMatrix{RGBf}`.

When `img` carries OME-Zarr pyramid levels, returns a `CompositePyramidSampler`
that selects the right resolution level and loads only the visible tile on each
zoom/pan event when wrapped in `Makie.Resampler`.  The `image!` overloads do
this automatically.

For images without a pyramid (or to force materialisation), returns a lazy
`MappedArray` / `colorview` backed by the original data.  Pass the result to
`collect` to materialise, or directly to `Makie.image` / `Makie.image!`.

Two modes selected automatically (or overridden via `colors`):

- **RGB mode** — channel labels exactly `["R","G","B"]`: lazy `colorview`, zero allocation.
- **Fluorescence mode** — all other configs: subsampled quantile pass to find
  per-channel clip values, then lazy `mappedarray` blend.

```julia
# Pyramid-aware (zoom-responsive) display
image!(ax, xen.images["morphology_focus"])  # handled automatically by image!

# Explicit composite for a cropped region
roi = view(xen, cx-500, cx+500, cy-500, cy+500)
rgb = composite(roi.images["morphology_focus"])
image!(ax, rgb)

# Explicit colours
rgb = composite(xen.images["morphology_focus"],
                colors=[:blue, :green, :red, :cyan])
```
"""
function composite(img::SpatialImage;
                   colors = nothing,
                   clip::Real = 0.999)
    ch_labels = channels(img)
    isempty(ch_labels) && error("composite: image has no channel axis")
    n = length(ch_labels)

    is_rgb = colors === nothing &&
             map(lowercase, ch_labels) == ["r", "g", "b"]
    resolved = _resolve_composite_colors(colors, n)

    if !isempty(img.pyramid)
        return CompositePyramidSampler(img; colors=colors, clip=clip)
    end

    return is_rgb ? _composite_rgb(img.data) :
                    _composite_fluor(img.data, resolved, Float64(clip))
end

# SpatialElementView overload — crop first, then lazy-blend the region.
# Always uses the non-pyramid path since we've already bounded the region.
function composite(v::SpatialElementView{<:SpatialImage};
                   colors = nothing,
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

    ch_labels = channels(img)
    n = length(ch_labels)
    is_rgb = colors === nothing && map(lowercase, ch_labels) == ["r", "g", "b"]
    resolved = _resolve_composite_colors(colors, n)

    return is_rgb ? _composite_rgb(cropped) : _composite_fluor(cropped, resolved, Float64(clip))
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
function _composite_fluor(data::AbstractArray{<:Real, 3},
                           colors::Vector{RGBf},
                           clip::Float64)
    n_ch = size(data, 1)

    his = Float32[_subsample_quantile(view(data, i, :, :), clip) for i in 1:n_ch]

    normed = [mappedarray(x -> clamp(Float32(x) / his[i], 0f0, 1f0),
                          view(data, i, :, :))
              for i in 1:n_ch]

    cr = Float32[c.r for c in colors]
    cg = Float32[c.g for c in colors]
    cb = Float32[c.b for c in colors]

    blend = let cr = cr, cg = cg, cb = cb, n = n_ch
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
function _pyramid_image!(ax, sampler::CompositePyramidSampler; kwargs...)
    nx, ny = size(sampler)   # (nx, ny) finest-resolution

    # Background: fixed-extent low-res overview — never updated.
    # Anchors the data bounds so reset_limits! / Ctrl+click always resets to
    # (1..nx) × (1..ny) regardless of where the detail layer is focused.
    ov_res  = min(nx, ny, 128)
    ov_tile = sampler(LinRange(1f0, Float32(nx), ov_res),
                      LinRange(1f0, Float32(ny), ov_res))
    bg = Makie.image!(ax, (1f0, Float32(nx)), (1f0, Float32(ny)), ov_tile; kwargs...)
    translate!(bg, 0, 0, -1)   # behind the detail layer

    # Detail: starts as a higher-res overview; updated to the visible tile on release.
    init_res  = min(nx, ny, 512)
    init_tile = sampler(LinRange(1f0, Float32(nx), init_res),
                        LinRange(1f0, Float32(ny), init_res))
    detail = Makie.image!(ax, (1f0, Float32(nx)), (1f0, Float32(ny)), init_tile; kwargs...)

    # Mirror HeatmapShader's debounce: advance slow_limits only when no mouse
    # buttons are held — tile loads on release, not during drag.
    events      = ax.scene.events
    slow_limits = Observable(ax.finallimits[])

    onany(ax.finallimits, events.mousebutton, events.keyboardbutton) do lims, _, _
        isempty(events.mousebuttonstate) || return
        last = slow_limits[]
        minimum(lims) ≈ minimum(last) && widths(lims) ≈ widths(last) && return
        slow_limits[] = lims
    end

    on(slow_limits) do lims
        xmin, ymin = Float32.(minimum(lims))
        xmax, ymax = Float32.(maximum(lims))

        x1 = clamp(xmin, 1f0, Float32(nx));  x2 = clamp(xmax, 1f0, Float32(nx))
        y1 = clamp(ymin, 1f0, Float32(ny));  y2 = clamp(ymax, 1f0, Float32(ny))
        (x2 <= x1 || y2 <= y1) && return

        res_x = clamp(round(Int, x2 - x1), 64, 2048)
        res_y = clamp(round(Int, y2 - y1), 64, 2048)
        tile  = sampler(LinRange(x1, x2, res_x), LinRange(y1, y2, res_y))

        detail[1][] = Makie.EndPoints{Float32}(x1, x2)
        detail[2][] = Makie.EndPoints{Float32}(y1, y2)
        detail[3][] = tile
    end

    return detail
end

"""
    image(img::SpatialImage; colors=nothing, clip=0.999, kwargs...)
    image!(ax, img::SpatialImage; colors=nothing, clip=0.999, kwargs...)

Render a `SpatialImage` as an RGB composite.

When pyramid levels are present (OME-Zarr data), uses `_pyramid_image!` which
subscribes to `ax.finallimits` for zoom-responsive tile loading — only the
visible region at the appropriate resolution is read from disk on each
zoom/pan event.

For images without a pyramid, calls `composite(img)` and collects at render time.
"""
function Makie.image(img::SpatialImage;
                     colors = nothing, clip = 0.999,
                     axis = (;), figure = (;), kwargs...)
    if !isempty(img.pyramid)
        fig = Figure(; figure...)
        ax  = Axis(fig[1, 1]; yreversed=true, axis...)
        sampler = CompositePyramidSampler(img; colors=colors, clip=clip)
        plt = _pyramid_image!(ax, sampler; kwargs...)
        return Makie.FigureAxisPlot(fig, ax, plt)
    else
        rgb = composite(img; colors=colors, clip=clip)
        ny, nx = size(rgb)
        return Makie.image(1..nx, 1..ny, permutedims(collect(rgb));
                           axis=axis, figure=figure, kwargs...)
    end
end

function Makie.image!(ax, img::SpatialImage; colors=nothing, clip=0.999, kwargs...)
    if !isempty(img.pyramid)
        sampler = CompositePyramidSampler(img; colors=colors, clip=clip)
        return _pyramid_image!(ax, sampler; kwargs...)
    else
        rgb = composite(img; colors=colors, clip=clip)
        ny, nx = size(rgb)
        return Makie.image!(ax, 1..nx, 1..ny, permutedims(collect(rgb)); kwargs...)
    end
end

function Makie.image!(ax, v::SpatialElementView{<:SpatialImage}; colors=nothing, clip=0.999, kwargs...)
    rgb = composite(v; colors=colors, clip=clip)
    e   = v.extent
    ny, nx = size(rgb)
    result = Makie.image!(ax, e.xmin..e.xmax, e.ymin..e.ymax,
                          permutedims(collect(rgb)); kwargs...)
    limits!(ax, e.xmin, e.xmax, e.ymin, e.ymax)
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
        heatmap!(ax, img; channel=channel, colormap=image_colormap)
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
