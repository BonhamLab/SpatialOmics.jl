# SpatialViz/src/recipes/labels.jl
#
# convert_arguments for SpatialLabels → heatmap
#
# SpatialLabels are integer masks (0 = background, positive = instance ID).
# We cast to Float32 so Makie's colormap machinery accepts them; the caller
# can set `colorrange`, `colormap`, and `lowclip = :transparent` as needed.
#
# Usage:
#   heatmap(lbl)
#   heatmap!(ax, lbl; colormap=:tab20, colorrange=(1, max_label))

"""
    Makie.convert_arguments(P, lbl::SpatialLabels)

Convert a `SpatialLabels` to heatmap arguments. The label array is cast to
`Float32` and transposed to Makie's (nx, ny) convention. For OME-Zarr label
arrays stored as (1, ny, nx), the singleton channel axis is dropped.
"""
function Makie.convert_arguments(P::Type{<:Heatmap}, lbl::SpatialLabels)
    data = lbl.data
    # Drop singleton leading dimension (e.g. (1, ny, nx) → (ny, nx))
    arr_2d = if ndims(data) == 3 && size(data, 1) == 1
        view(data, 1, :, :)
    elseif ndims(data) == 2
        data
    else
        error("SpatialLabels: expected 2-D or (1,ny,nx) array, got size $(size(data))")
    end

    ny, nx = size(arr_2d)
    mat = permutedims(Float32.(collect(arr_2d)))  # (nx, ny) for Makie
    return Makie.convert_arguments(P, 1:nx, 1:ny, mat)
end
