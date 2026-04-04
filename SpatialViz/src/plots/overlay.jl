# SpatialViz/src/plots/overlay.jl
# Composite multi-layer figures (image + points, image + labels, etc.)

"""
    overlay_points_on_image(img::SpatialImage, pts::SpatialPoints;
                             color = :gene,
                             markersize::Real = 2.0,
                             image_colormap = :grays,
                             backend::Type{<:VisualizationBackend} = MakieBackend,
                             kwargs...)

Render `pts` scatter overlaid on `img` in a single figure.
The image is shown in `image_colormap`; points are coloured by `color`.
"""
function overlay_points_on_image(
    img::SpatialImage,
    pts::SpatialPoints;
    color = :gene,
    markersize::Real = 2.0,
    image_colormap = :grays,
    backend::Type{<:VisualizationBackend} = MakieBackend,
    kwargs...,
)
    # TODO Phase 1: implement for MakieBackend
    error("overlay_points_on_image: not yet implemented (Phase 1 deliverable)")
end
