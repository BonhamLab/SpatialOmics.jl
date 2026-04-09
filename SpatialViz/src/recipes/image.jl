# SpatialViz/src/recipes/image.jl
#
# convert_arguments for SpatialImage → heatmap
#
# Wraps multi-resolution images in ImagePyramidSampler + Makie.Resampler for
# zoom-responsive display. Single-resolution images fall back to a plain
# heatmap.
#
# Usage (after loading a backend):
#   heatmap(img)              # channel 1, grays colormap
#   heatmap(img; channel=2)   # NOT supported via convert_arguments —
#                             #   use the helper below instead:
#   spatial_image(img, channel=2)   # returns (x_range, y_range, Resampler)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""
    _channel_view(img::SpatialImage, channel::Int=1) -> AbstractArray{T,2}

Return a lazy 2-D view of `channel` from an OME-Zarr (c, y, x) image array.
For 2-D inputs the array is returned as-is.
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
# convert_arguments — SpatialImage → heatmap (channel 1)
# ---------------------------------------------------------------------------

"""
    Makie.convert_arguments(P, img::SpatialImage)

Converts a `SpatialImage` to heatmap arguments. If the image carries pyramid
levels they are wrapped in `Makie.Resampler` for zoom-responsive lazy loading.
Channel 1 is used; call `spatial_image(img; channel=k)` to select a different
channel.
"""
function Makie.convert_arguments(P::Type{<:Heatmap}, img::SpatialImage)
    return _image_to_heatmap_args(P, img, 1)
end

function _image_to_heatmap_args(P, img::SpatialImage, channel::Int)
    data_2d = _channel_view(img.data, channel)
    ny, nx   = size(data_2d)

    if !isempty(img.pyramid)
        sampler   = ImagePyramidSampler(img, channel)
        resampler = Makie.Resampler(sampler; lowres_background=true)
        # Pass explicit spatial extents so axis limits match pixel coordinates.
        return Makie.convert_arguments(P, 1:nx, 1:ny, resampler)
    else
        # Single-resolution: materialise and transpose to (nx, ny) for heatmap.
        # heatmap(xs, ys, data) expects data[i,j] at (xs[i], ys[j]).
        mat = collect(Float32.(data_2d))   # materialise DiskArray once
        return Makie.convert_arguments(P, 1:nx, 1:ny, permutedims(mat))
    end
end

# ---------------------------------------------------------------------------
# heatmap overloads — SpatialElementView{SpatialImage} (view with extent)
#
# Renders the full image (Resampler selects the right pyramid level) then
# constrains the axis to the view extent, so the caller sees only that region.
# ---------------------------------------------------------------------------

function Makie.heatmap!(ax, v::SpatialElementView{<:SpatialImage}; channel::Int=1, kwargs...)
    result = heatmap!(ax, v.parent; channel=channel, kwargs...)
    e = v.extent
    limits!(ax, e.xmin, e.xmax, e.ymin, e.ymax)
    return result
end

function Makie.heatmap(v::SpatialElementView{<:SpatialImage}; channel::Int=1, kwargs...)
    result = heatmap(v.parent; channel=channel, kwargs...)
    e = v.extent
    limits!(result[2], e.xmin, e.xmax, e.ymin, e.ymax)
    return result
end

# ---------------------------------------------------------------------------
# heatmap overloads — dispatch on SpatialImage with channel kwarg
# ---------------------------------------------------------------------------

"""
    heatmap(img::SpatialImage; channel::Int=1, kwargs...)
    heatmap!(ax, img::SpatialImage; channel::Int=1, kwargs...)

Plot a `SpatialImage` as a heatmap, selecting `channel` from the
`(c, y, x)` array. When pyramid levels are present, `Makie.Resampler`
is used for zoom-responsive lazy tile loading.

```julia
fig, ax, plt = heatmap(img; channel=1, colormap=:grays)
heatmap!(ax, img; channel=2, colormap=:viridis)
```
"""
function Makie.heatmap(img::SpatialImage; channel::Int=1, kwargs...)
    args = _image_to_heatmap_args(Heatmap, img, channel)
    return heatmap(args...; kwargs...)
end

function Makie.heatmap!(ax, img::SpatialImage; channel::Int=1, kwargs...)
    args = _image_to_heatmap_args(Heatmap, img, channel)
    return heatmap!(ax, args...; kwargs...)
end
