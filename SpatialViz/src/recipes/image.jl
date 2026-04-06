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
    ny, nx = size(s)
    return ((1f0, Float32(nx)), (1f0, Float32(ny)), s)
end

function Makie.convert_arguments(P::Type{<:Heatmap}, x, y, s::ImagePyramidSampler)
    ny, nx = size(s)
    return ((1f0, Float32(nx)), (1f0, Float32(ny)), s)
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
# spatial_image — thin helper that respects channel kwarg
# ---------------------------------------------------------------------------

"""
    spatial_image(img::SpatialImage; channel::Int=1, kwargs...)

Return a `(Figure, Axis, plot)` tuple showing `img` in a new figure.
This thin wrapper exists only to allow `channel` selection; for compositing
with other layers use the `convert_arguments` dispatch directly.

```julia
fig, ax, plt = spatial_image(ds.images["morphology_focus"]; channel=1)
```
"""
function spatial_image(img::SpatialImage; channel::Int=1, kwargs...)
    args = _image_to_heatmap_args(Heatmap, img, channel)
    return heatmap(args...; kwargs...)
end

function spatial_image!(ax, img::SpatialImage; channel::Int=1, kwargs...)
    args = _image_to_heatmap_args(Heatmap, img, channel)
    return heatmap!(ax, args...; kwargs...)
end
