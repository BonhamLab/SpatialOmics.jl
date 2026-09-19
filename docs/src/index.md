# SpatialOmics.jl

A Julia library for loading, representing, and analysing spatial transcriptomics
data. It provides a common data model for multi-modal spatial experiments —
transcripts, cell boundaries, tissue images, segmentation masks, and
expression matrices — alongside lazy spatial views, a multi-FOV coordinate
system graph, explicit persistence, and import support for SpatialData Zarr
stores.

## Installation

```julia
using Pkg
Pkg.add("SpatialOmics")
```

## Quick start

```julia
using CairoMakie   # load a Makie backend before plotting
using SpatialOmics

# Load a native store or a supported SpatialData Zarr store
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
image!(ax,   scaleminmax(channel(images(roi, "morphology_focus"), 1)))
poly!(ax,    shapes(roi, "cell_boundaries");  color=:transparent, strokecolor=:cyan)
scatter!(ax, points(roi, "transcripts");      markersize=1, color=(:red, 0.3))
tightlimits!(ax)
fig
```

## Navigation

- **[Explanation](@ref "The data model")** — Why things are designed the way they
  are: the data model, coordinate system graph, and lazy view semantics.
- **[Tutorials](@ref "Tutorials")** — Executable core lessons plus
  pre-rendered workflows using public technology datasets.
- **[Reference](@ref "Dataset")** — Complete API documentation for all exported
  functions and types.
- **[Guides](@ref "Quickstart")** — Task-oriented walkthroughs: loading data,
  building plots, working with CosMx exports.
