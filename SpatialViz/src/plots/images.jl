# SpatialViz/src/plots/images.jl
# Raster image layer rendering (DAPI, H&E, IF channels).

"""
    plot_image(img::SpatialImage;
               channel::Union{Int,Nothing} = nothing,
               colormap = :grays,
               contrast_limits::Union{Nothing, Tuple{Real,Real}} = nothing,
               backend::Type{<:VisualizationBackend} = MakieBackend,
               kwargs...)

Render a single channel of `img` as a 2-D heatmap/image. If `channel` is
`nothing` and `img` has multiple channels, the first channel is used.
`contrast_limits` pins the colormap range; pass `nothing` for auto-scaling.
"""
function plot_image(
    img::SpatialImage;
    channel::Union{Int,Nothing} = nothing,
    colormap = :grays,
    contrast_limits::Union{Nothing, Tuple{Real,Real}} = nothing,
    backend::Type{<:VisualizationBackend} = MakieBackend,
    kwargs...,
)
    # TODO Phase 1: implement for MakieBackend
    error("plot_image: not yet implemented (Phase 1 deliverable)")
end
