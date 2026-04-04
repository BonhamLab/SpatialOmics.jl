# SpatialViz/src/spatial_plot.jl
# Top-level spatial_plot() dispatcher — the primary public API entry-point.

"""
    spatial_plot(data::SpatialDataset;
                 backend::Type{<:VisualizationBackend} = MakieBackend,
                 interactive::Bool = true,
                 layers::Vector{String} = ["points", "images"],
                 color::Union{Symbol,String,Nothing} = nothing,
                 colormap = :viridis,
                 kwargs...)

Render a spatial visualization of `data` using the specified `backend`.

Layer names in `layers` correspond to keys in the relevant element dictionaries
of `data` (e.g. `"morphology"` for `data.images["morphology"]`), or the
special strings `"images"`, `"points"`, `"labels"`, `"shapes"` which render
all elements of that type.

`color` may be a column name in the observation table, a gene name, or `nothing`
to use a default coloring.

# Examples
```julia
# Simple high-level call — renders transcript scatter overlaid on DAPI image
spatial_plot(ds)

# Explicit backend and layer control
spatial_plot(ds;
    backend = WGLMakieBackend,
    interactive = true,
    layers = ["morphology", "transcripts"],
    color = "cell_type",
    colormap = :tab20)
```
"""
function spatial_plot(
    data::SpatialDataset;
    backend::Type{<:VisualizationBackend} = MakieBackend,
    interactive::Bool = true,
    layers::Vector{String} = ["points", "images"],
    color::Union{Symbol,String,Nothing} = nothing,
    colormap = :viridis,
    kwargs...,
)
    renderer = create_renderer(backend, interactive)
    return render_spatial(renderer, data, layers; color = color, colormap = colormap, kwargs...)
end
