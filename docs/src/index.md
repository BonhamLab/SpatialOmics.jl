# SpatialOmics.jl

A Julia library for loading, representing, and analysing spatial transcriptomics
data. It provides a common data model for multi-modal spatial experiments —
transcripts, cell boundaries, morphology images, segmentation masks, and
expression matrices — alongside lazy spatial views, a multi-FOV coordinate
system graph, and SpatialData OME-Zarr interoperability with Python tools.

## Installation

```julia
using Pkg
Pkg.add("SpatialOmics")
```

## Quick start

```julia
using CairoMakie   # load a Makie backend before plotting
using SpatialOmics

# Load from SpatialData OME-Zarr (Xenium, CosMx, Visium, …)
ds = read(SpatialDataZarr(), "/path/to/experiment.zarr")

# Inspect structure
keys(elements(ds))
top_features(points(ds, "transcripts"), 10)

# Define a region of interest — lazy, no data copied
ext = SpatialExtent(4000.0, 5000.0, 1000.0, 2000.0; coord_system="global")
roi = view(ds, ext)

# Build a composite panel using standard Makie verbs
fig = Figure(size=(600, 600))
ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true)
heatmap!(ax, images(roi, "morphology_focus"); channel=1, colormap=:grays)
poly!(ax,    shapes(roi, "cell_boundaries");  color=:transparent, strokecolor=:cyan)
scatter!(ax, points(roi, "transcripts");      markersize=1, color=(:red, 0.3))
tightlimits!(ax)
fig
```

## Navigation

- **[Explanation](@ref "The data model")** — Why things are designed the way they
  are: the data model, coordinate system graph, and lazy view semantics.
- **[Reference](@ref "Dataset")** — Complete API documentation for all exported
  functions and types.
- **[Guides](@ref "Quickstart")** — Task-oriented walkthroughs: loading data,
  building plots, working with CosMx exports.
