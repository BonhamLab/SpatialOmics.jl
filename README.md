# SpatialOmics.jl

[![CI](https://github.com/BonhamLab/SpatialOmics.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/BonhamLab/SpatialOmics.jl/actions/workflows/ci.yml)
<!-- [![Docs (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://BonhamLab.github.io/SpatialOmics.jl/stable/) -->

[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://BonhamLab.github.io/SpatialOmics.jl/dev/)
[![Project Status: WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

A Julia library for loading, representing, and analysing spatial transcriptomics data.
It provides a common data model for multi-modal spatial experiments — transcripts, cell boundaries, tissue images, segmentation masks, and expression matrices — alongside lazy spatial views, a multi-FOV coordinate system graph, explicit persistence, and import support for SpatialData Zarr stores.

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

See the [documentation](https://BonhamLab.github.io/SpatialOmics.jl/dev/) for a full API reference, explanations of the data model and coordinate system graph, and platform-specific guides.
