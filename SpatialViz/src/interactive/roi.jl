# SpatialViz/src/interactive/roi.jl
# Region-of-interest (ROI) selection tools.

"""
    roi_select(ds::SpatialDataset;
               tool::Symbol = :polygon,
               backend::Type{<:VisualizationBackend} = MakieBackend)
        -> SpatialDataset

Launch an interactive ROI selection overlay on the current spatial plot.
Returns a new `SpatialDataset` containing only the observations inside the
selected region.

`tool` options (planned): `:polygon`, `:lasso`, `:rectangle`, `:circle`.
"""
function roi_select(
    ds::SpatialDataset;
    tool::Symbol = :polygon,
    backend::Type{<:VisualizationBackend} = MakieBackend,
)
    error("roi_select: not yet implemented (Phase 2 deliverable)")
end

# Convenience aliases
lasso_select(ds::SpatialDataset; kwargs...) = roi_select(ds; tool = :lasso, kwargs...)
rectangle_select(ds::SpatialDataset; kwargs...) = roi_select(ds; tool = :rectangle, kwargs...)
