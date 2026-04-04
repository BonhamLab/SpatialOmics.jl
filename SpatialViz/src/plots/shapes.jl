# SpatialViz/src/plots/shapes.jl
# Vector shape layer rendering (cell boundary polygons, tissue region outlines).

"""
    plot_shapes(shp::SpatialShapes;
                color = :cell_type,
                colormap = :tab20,
                linewidth::Real = 0.5,
                fill_alpha::Real = 0.0,
                backend::Type{<:VisualizationBackend} = MakieBackend,
                kwargs...)

Render cell boundary polygons or tissue regions. By default outlines only
(`fill_alpha=0`); set `fill_alpha > 0` for filled polygons.
"""
function plot_shapes(
    shp::SpatialShapes;
    color = :cell_type,
    colormap = :tab20,
    linewidth::Real = 0.5,
    fill_alpha::Real = 0.0,
    backend::Type{<:VisualizationBackend} = MakieBackend,
    kwargs...,
)
    # TODO Phase 2: implement
    error("plot_shapes: not yet implemented (Phase 2 deliverable)")
end
