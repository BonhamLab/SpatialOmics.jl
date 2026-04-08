"""
    SpatialViz.jl

Makie-based visualization for the SpatialOmics ecosystem.

Extends Makie's `convert_arguments` to dispatch directly on SpatialOmicsBase
element types — users plot spatial data with the same Makie verbs they already
know (`heatmap`, `scatter`, `poly`, etc.).

**Backend choice is left to the user.** Load a Makie backend *before* using
SpatialViz:
```julia
using GLMakie      # native window, recommended for interactive exploration
# or
using CairoMakie   # static output, recommended for notebooks / publication
# or
using WGLMakie     # browser-based, for Pluto / Jupyter

using SpatialViz
```

SpatialViz depends only on the abstract `Makie` package so it is
backend-agnostic and will not pull in GL/WGL/Cairo deps itself.

## Pyramid-aware image display
`SpatialImage` objects that carry multi-resolution pyramid levels (loaded from
OME-Zarr stores) are wrapped in `ImagePyramidSampler` and passed to
`Makie.Resampler`, giving zoom-responsive lazy tile loading at no extra cost.

## Coordinate-based selection
`SpatialExtent` and `SpatialView` (from SpatialOmicsBase) can be used to crop
elements before plotting — see the `crop` function.
"""
module SpatialViz

using Makie
using GeometryBasics
using SpatialOmicsBase
using ImageCore: colorview, RGB
using MappedArrays: mappedarray

# ---------------------------------------------------------------------------
# Re-export key types users will need
# ---------------------------------------------------------------------------
export SpatialExtent, SpatialElementView, SpatialDatasetView
export extent, intersects, crop
export ImagePyramidSampler

# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
export spatial_panel, composite

# ---------------------------------------------------------------------------
# Includes — one file per element type, plus composite helpers
# ---------------------------------------------------------------------------
include("recipes/image.jl")
include("recipes/points.jl")
include("recipes/shapes.jl")
include("recipes/labels.jl")
include("recipes/composite.jl")

end # module SpatialViz
