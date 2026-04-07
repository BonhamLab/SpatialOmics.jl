# SpatialViz/src/recipes/composite.jl
#
# Composite plot helpers — functions that create multi-layer figures by
# combining the primitive convert_arguments recipes above.
#
# These are thin conveniences; users can always compose layers manually via
# the standard Makie API (heatmap! + scatter! etc.).

"""
    spatial_panel(ds::SpatialDataset;
                  image_key="morphology_focus", channel=1,
                  points_key=nothing,
                  shapes_key=nothing,
                  labels_key=nothing,
                  image_colormap=:grays,
                  points_color=:red,
                  points_markersize=2,
                  shapes_color=:transparent,
                  shapes_strokecolor=:white,
                  shapes_strokewidth=0.5,
                  view::Union{SpatialDatasetView,Nothing}=nothing,
                  kwargs...)

Create a `Figure` with a single `Axis` showing the requested layers from `ds`.

Layers are rendered in order: image → labels → shapes → points, so points
always appear on top. Pass `nothing` for any key to skip that layer.

Returns `(fig, ax)`.
"""
function spatial_panel(
    ds::SpatialDataset;
    image_key::Union{String,Nothing}  = _first_key(ds.images),
    channel::Int                      = 1,
    points_key::Union{String,Nothing} = _first_key(ds.points),
    shapes_key::Union{String,Nothing} = nothing,
    labels_key::Union{String,Nothing} = nothing,
    image_colormap                    = :grays,
    points_color                      = :crimson,
    points_markersize::Real           = 2,
    shapes_color                      = :transparent,
    shapes_strokecolor                = :white,
    shapes_strokewidth::Real          = 0.5,
    view::Union{SpatialDatasetView,Nothing} = nothing,
    figure_kwargs                     = (;),
    axis_kwargs                       = (;),
    kwargs...,
)
    fig = Figure(; figure_kwargs...)
    ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true, axis_kwargs...)

    # --- image layer ---
    if image_key !== nothing && haskey(ds.images, image_key)
        img = ds.images[image_key]
        heatmap!(ax, img; channel=channel, colormap=image_colormap)
    end

    # --- labels layer ---
    if labels_key !== nothing && haskey(ds.labels, labels_key)
        heatmap!(ax, ds.labels[labels_key]; colormap=:tab20)
    end

    # --- shapes layer ---
    if shapes_key !== nothing && haskey(ds.shapes, shapes_key)
        shp = view !== nothing ? view.shapes[shapes_key] : ds.shapes[shapes_key]
        poly!(ax, shp;
              color=shapes_color,
              strokecolor=shapes_strokecolor,
              strokewidth=shapes_strokewidth)
    end

    # --- points layer ---
    if points_key !== nothing && haskey(ds.points, points_key)
        pts = view !== nothing ? view.points[points_key] : ds.points[points_key]
        scatter!(ax, pts;
                 color=points_color,
                 markersize=points_markersize)
    end

    # Constrain axis to the view extent when one is provided; otherwise fit all layers.
    if view !== nothing
        e = view.extent
        limits!(ax, e.xmin, e.xmax, e.ymin, e.ymax)
    else
        tightlimits!(ax)
    end

    return fig, ax
end

_first_key(d::Dict) = isempty(d) ? nothing : first(keys(d))
