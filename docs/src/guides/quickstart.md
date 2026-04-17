# Quickstart

## Loading data

```julia
using SpatialOmics
import SpatialOmics as SO

# SpatialData OME-ZARR (Xenium, CosMx, Visium, MERFISH)
ds = SO.read(SO.Zarr(), "/path/to/experiment.zarr")

# Platform-specific loaders (un-exported; use SO.* prefix)
xen = SO.load_xenium("/path/to/xenium_output/")
cos = SO.load(SO.CosMxReader(), "/path/to/cosmx_export/")

# Round-trip back to SpatialData OME-ZARR
SO.write(ds, "/path/to/output.zarr", SO.Zarr())
```

## Inspecting a dataset

```julia
# Element dictionaries
keys(ds.images)   # e.g. ["morphology_focus"]
keys(ds.points)   # e.g. ["transcripts"]
keys(ds.shapes)   # e.g. ["cell_boundaries"]
keys(ds.tables)   # e.g. ["cell_by_gene"]

# Element access
img = ds.images["morphology_focus"]
pts = ds.points["transcripts"]

# Bounding box of any element
ext = extent(pts)   # SpatialExtent{Float32}
```

## Defining a region of interest

`SpatialExtent` is a typed bounding box. `view(ds, ext)` returns a
`SpatialDatasetView` — a lazy pointer into the dataset that applies the
spatial filter only when data are actually needed:

```julia
# Centre on the transcript cloud
cx = (ext.xmin + ext.xmax) / 2
cy = (ext.ymin + ext.ymax) / 2
roi = view(ds, SpatialExtent(cx-500, cx+500, cy-500, cy+500))

# Element access returns SpatialElementView — still no data copied
pts_view = roi.points["transcripts"]
```

## Visualization

Load a Makie backend **before** `SpatialViz`. The package is backend-agnostic.

```julia
using GLMakie     # interactive; or CairoMakie (static), WGLMakie (browser)
using SpatialViz
```

Standard Makie verbs are extended to accept spatial types directly:

```julia
# Full slide — Makie.Resampler selects the correct pyramid level on each zoom
heatmap(ds.images["morphology_focus"]; channel=1, colormap=:grays)

# With a region view — axis is constrained to the ROI automatically
fig = Figure(size=(600, 600))
ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true)
heatmap!(ax, roi.images["morphology_focus"]; channel=1, colormap=:grays)
poly!(ax,    roi.shapes["cell_boundaries"];  color=:transparent, strokecolor=:cyan)
scatter!(ax, roi.points["transcripts"];      markersize=1, color=(:red, 0.3))
fig
```

See the [Visualization guide](@ref) for details on pyramid levels, channel
selection, and multi-panel layouts.
