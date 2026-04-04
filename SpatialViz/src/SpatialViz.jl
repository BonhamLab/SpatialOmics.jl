"""
    SpatialViz.jl

Visualization framework for the STX_DEV spatial transcriptomics ecosystem.
Implements a multi-backend architecture dispatching on subtypes of
`VisualizationBackend`: `MakieBackend` for native Julia static/interactive
plots, `WGLMakieBackend` / `BonitoBackend` for browser-based deployment, and
`NapariBackend` for Python Napari via PythonCall.jl.

The single public entry-point `spatial_plot()` handles all rendering; layer
types, color maps, and interactivity options are controlled via keyword args.
"""
module SpatialViz

using Makie
using WGLMakie
using Bonito
using PythonCall
using SpatialOmicsBase

# ---------------------------------------------------------------------------
# Exports — backend types
# ---------------------------------------------------------------------------

export VisualizationBackend
export MakieBackend, WGLMakieBackend, BonitoBackend, NapariBackend

# ---------------------------------------------------------------------------
# Exports — primary plotting API
# ---------------------------------------------------------------------------

export spatial_plot
export plot_points, plot_image, plot_labels, plot_shapes
export overlay_points_on_image

# ---------------------------------------------------------------------------
# Exports — interactive annotation tools
# ---------------------------------------------------------------------------

export annotation_tool
export roi_select, lasso_select, rectangle_select
export save_annotations, load_annotations

# ---------------------------------------------------------------------------
# Exports — color utilities
# ---------------------------------------------------------------------------

export categorical_palette, continuous_colormap, colorblind_safe_palette

# ---------------------------------------------------------------------------
# Includes
# ---------------------------------------------------------------------------

include("backends/backends.jl")
include("backends/makie_backend.jl")
include("backends/wglmakie_backend.jl")
include("backends/napari_backend.jl")
include("plots/points.jl")
include("plots/images.jl")
include("plots/labels.jl")
include("plots/shapes.jl")
include("plots/overlay.jl")
include("interactive/annotation.jl")
include("interactive/roi.jl")
include("colors/palettes.jl")
include("spatial_plot.jl")

end # module SpatialViz
