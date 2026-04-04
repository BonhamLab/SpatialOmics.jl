# SpatialViz/src/backends/makie_backend.jl
# Native Makie (GLMakie / CairoMakie) backend.

"""
    MakieBackend

Renders spatial data using Makie.jl. Supports both static (CairoMakie) and
interactive (GLMakie) rendering. Default backend for SpatialViz.
"""
struct MakieBackend <: VisualizationBackend end

function create_renderer(::Type{MakieBackend}, interactive::Bool)
    # TODO Phase 1: return a Makie Figure
    error("MakieBackend.create_renderer: not yet implemented (Phase 1 deliverable)")
end

function render_spatial(renderer, data::SpatialDataset, layers::Vector{String};
                        kwargs...)
    # TODO Phase 1: dispatch on layer type strings and call plot_* helpers
    error("MakieBackend.render_spatial: not yet implemented (Phase 1 deliverable)")
end
