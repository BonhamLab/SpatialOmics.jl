# SpatialViz/src/plots/labels.jl
# Integer label image rendering (cell/nucleus segmentation masks).

"""
    plot_labels(lbl::SpatialLabels;
                colormap = :tab20,
                alpha::Real = 0.5,
                show_boundaries::Bool = true,
                backend::Type{<:VisualizationBackend} = MakieBackend,
                kwargs...)

Render segmentation labels as a coloured overlay. When `show_boundaries=true`,
draw only the outlines of labelled regions rather than solid fills.
"""
function plot_labels(
    lbl::SpatialLabels;
    colormap = :tab20,
    alpha::Real = 0.5,
    show_boundaries::Bool = true,
    backend::Type{<:VisualizationBackend} = MakieBackend,
    kwargs...,
)
    # TODO Phase 2: implement
    error("plot_labels: not yet implemented (Phase 2 deliverable)")
end
