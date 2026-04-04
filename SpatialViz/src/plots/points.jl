# SpatialViz/src/plots/points.jl
# Scatter/point layer rendering for transcript and cell centroid data.

"""
    plot_points(pts::SpatialPoints;
                color = :gene,
                colormap = :tab20,
                markersize::Real = 2.0,
                alpha::Real = 0.8,
                backend::Type{<:VisualizationBackend} = MakieBackend,
                kwargs...)

Render `pts` as a 2-D scatter plot. `color` may be a column name in
`pts.features` (categorical or continuous) or a fixed color value.

Returns a backend-specific figure/axis object.
"""
function plot_points(
    pts::SpatialPoints;
    color = :gene,
    colormap = :tab20,
    markersize::Real = 2.0,
    alpha::Real = 0.8,
    backend::Type{<:VisualizationBackend} = MakieBackend,
    kwargs...,
)
    # TODO Phase 1: implement for MakieBackend
    error("plot_points: not yet implemented (Phase 1 deliverable)")
end
