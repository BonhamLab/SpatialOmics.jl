# SpatialViz/src/backends/backends.jl
# Abstract VisualizationBackend and the renderer factory.

"""
    VisualizationBackend

Abstract supertype for all visualization backends. Concrete subtypes implement
`create_renderer(backend, interactive)` and `render_spatial(renderer, data, layers; kwargs...)`.
"""
abstract type VisualizationBackend end

"""
    create_renderer(backend::Type{<:VisualizationBackend}, interactive::Bool)

Return a backend-specific renderer object (Figure, Canvas, etc.).
Must be implemented by each concrete `VisualizationBackend` subtype.
"""
function create_renderer end

"""
    render_spatial(renderer, data::SpatialDataset, layers::Vector{String}; kwargs...)

Draw the requested `layers` of `data` onto `renderer`.
Must be implemented by each concrete `VisualizationBackend` subtype.
"""
function render_spatial end
