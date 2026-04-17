# SpatialViz/src/recipes/image.jl
#
# heatmap overloads for SpatialImage (and SpatialElementView{SpatialImage}).
#
# Both dispatch through the same methods via two private accessors:
#   _img_parent(x)  → the underlying SpatialImage
#   _img_extent(x)  → SpatialExtent or nothing
#
# After plotting, limits! is applied when an extent is present.

# ---------------------------------------------------------------------------
# Private accessors — unify SpatialImage and SpatialElementView{SpatialImage}
# ---------------------------------------------------------------------------

_img_parent(img::SpatialImage)                        = img
_img_parent(v::SpatialElementView{<:SpatialImage})    = v.parent

_img_extent(::SpatialImage)                           = nothing
_img_extent(v::SpatialElementView{<:SpatialImage})    = v.extent

const SpatialImageLike = Union{SpatialImage, SpatialElementView{<:SpatialImage}}

# Apply axis limits when the input carries a spatial extent.
_apply_extent!(ax, ::Nothing) = nothing
_apply_extent!(ax, e::SpatialExtent) = limits!(ax, e.xmin, e.xmax, e.ymin, e.ymax)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""
    _channel_view(arr, channel=1) -> AbstractArray{T,2}

Return a lazy 2-D view of `channel` from a (c, y, x) array. 2-D arrays are
returned as-is.
"""
function _channel_view(arr::AbstractArray, channel::Int=1)
    ndims(arr) == 3 && return view(arr, channel, :, :)
    ndims(arr) == 2 && return arr
    error("_channel_view: expected 2-D or 3-D array, got $(ndims(arr))-D")
end

# ---------------------------------------------------------------------------
# convert_arguments for ImagePyramidSampler
#
# Makie's Resampler path calls:
#   convert_arguments(Heatmap, x, y, resampler)
#     → convert_arguments(Heatmap, x, y, resampler.data)   ← hits our override
#     → (EndPoints{Float32}(x...), EndPoints{Float32}(y...), resampler)
#
# Without this override Makie falls through to the generic AbstractMatrix path,
# which tries to materialise the full 17k×51k sampler → OOM.
# ---------------------------------------------------------------------------

function Makie.convert_arguments(P::Type{<:Heatmap}, s::ImagePyramidSampler)
    nx, ny = size(s)
    x1, x2 = s.x_range !== nothing ? s.x_range : (1.0, Float64(nx))
    y1, y2 = s.y_range !== nothing ? s.y_range : (1.0, Float64(ny))
    return ((Float32(x1), Float32(x2)), (Float32(y1), Float32(y2)), s)
end

function Makie.convert_arguments(P::Type{<:Heatmap}, x, y, s::ImagePyramidSampler)
    nx, ny = size(s)
    x1, x2 = s.x_range !== nothing ? s.x_range : (1.0, Float64(nx))
    y1, y2 = s.y_range !== nothing ? s.y_range : (1.0, Float64(ny))
    return ((Float32(x1), Float32(x2)), (Float32(y1), Float32(y2)), s)
end

# ---------------------------------------------------------------------------
# convert_arguments — SpatialImage → heatmap (channel 1, no kwarg path)
# ---------------------------------------------------------------------------

function Makie.convert_arguments(P::Type{<:Heatmap}, img::SpatialImage)
    return _image_to_heatmap_args(P, img, 1)
end

function _image_to_heatmap_args(P, img::SpatialImage, channel::Int)
    data_2d = _channel_view(img.data, channel)
    ny, nx   = size(data_2d)

    if !isempty(img.pyramid)
        sampler   = ImagePyramidSampler(img, channel)
        resampler = Makie.Resampler(sampler; lowres_background=true)
        nx_s, ny_s = size(sampler)
        x_ext = sampler.x_range !== nothing ? sampler.x_range : (1.0, Float64(nx_s))
        y_ext = sampler.y_range !== nothing ? sampler.y_range : (1.0, Float64(ny_s))
        return Makie.convert_arguments(P, x_ext, y_ext, resampler)
    else
        # Downsample to at most MAX_PX on the long axis if no pyramid is present.
        # Prevents materialising multi-GB arrays (e.g. full-slide CosMx morphology).
        MAX_PX = 4096
        if ny > MAX_PX || nx > MAX_PX
            ratio  = ny / nx
            ny_s   = ratio >= 1.0 ? min(ny, MAX_PX) : max(1, round(Int, min(nx, MAX_PX) * ratio))
            nx_s   = ratio >= 1.0 ? max(1, round(Int, ny_s / ratio))  : min(nx, MAX_PX)
            y_step = max(1, ny ÷ ny_s)
            x_step = max(1, nx ÷ nx_s)
            mat = collect(Float32.(data_2d[1:y_step:end, 1:x_step:end]))
        else
            mat = collect(Float32.(data_2d))
        end
        xr    = get(img.metadata, "x_range", nothing)
        yr    = get(img.metadata, "y_range", nothing)
        x_ext = xr !== nothing ? (Float64(xr[1]), Float64(xr[2])) : (1.0, Float64(nx))
        y_ext = yr !== nothing ? (Float64(yr[1]), Float64(yr[2])) : (1.0, Float64(ny))
        return Makie.convert_arguments(P, x_ext, y_ext, permutedims(mat))
    end
end

# ---------------------------------------------------------------------------
# heatmap / heatmap! — SpatialImage and SpatialElementView{SpatialImage}
# ---------------------------------------------------------------------------

"""
    heatmap(img; channel=1, kwargs...)
    heatmap!(ax, img; channel=1, kwargs...)

Plot a `SpatialImage` (or a `SpatialElementView{SpatialImage}`) as a heatmap,
selecting `channel` from the `(c, y, x)` array. When pyramid levels are
present, `Makie.Resampler` is used for zoom-responsive lazy tile loading.
When called on a view, the axis is constrained to the view's extent.

```julia
fig, ax, plt = heatmap(img; channel=1, colormap=:grays)
heatmap!(ax, img; channel=2, colormap=:viridis)
heatmap!(ax, roi.images["morphology_focus"]; channel=1, colormap=:grays)
```
"""
function Makie.heatmap(img::SpatialImageLike; channel::Int=1, kwargs...)
    args = _image_to_heatmap_args(Heatmap, _img_parent(img), channel)
    fig, ax, plt = heatmap(args...; kwargs...)
    _apply_extent!(ax, _img_extent(img))
    return Makie.FigureAxisPlot(fig, ax, plt)
end

function Makie.heatmap!(ax, img::SpatialImageLike; channel::Int=1, kwargs...)
    args = _image_to_heatmap_args(Heatmap, _img_parent(img), channel)
    result = heatmap!(ax, args...; kwargs...)
    _apply_extent!(ax, _img_extent(img))
    return result
end
