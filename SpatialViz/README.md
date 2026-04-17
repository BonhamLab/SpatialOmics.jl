# SpatialViz.jl

SpatialViz.jl is the visualization layer of the SpatialOmics.jl monorepo. It extends Makie.jl with spatial-omics-aware recipes for all `SpatialDataset` element types:

- **`heatmap` / `heatmap!`** — renders `SpatialImage` with automatic coordinate registration from `x_range`/`y_range` metadata. When OME-Zarr pyramid levels are present, uses `Makie.Resampler` + `ImagePyramidSampler` for zoom-responsive lazy tile loading; falls back to downsampled materialize for non-pyramidal arrays.
- **`scatter` / `scatter!`** — renders `SpatialPoints` with optional feature-column colouring.
- **`poly` / `poly!`** — renders `SpatialShapes` (polygons, circles) with optional per-cell colouring.

All recipes accept both the underlying element type (`SpatialImage`, `SpatialPoints`, `SpatialShapes`) and lazy `SpatialElementView{T}` objects returned by `view(ds, SpatialExtent(...))`, automatically applying axis limits to the view's extent.

SpatialViz depends only on the abstract `Makie` package — users load a concrete backend (`CairoMakie`, `GLMakie`, or `WGLMakie`) before `using SpatialViz`.
