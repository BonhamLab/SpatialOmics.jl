# SpatialOmics.jl

A modular, high-performance Julia ecosystem for spatial transcriptomics analysis.

## Packages

| Package | Role |
|---------|------|
| [`SpatialOmicsBase`](@ref) | Core types — `SpatialDataset`, elements, coordinate systems, lazy views |
| [`SpatialIO`](@ref) | All I/O: platform readers, SpatialData interop, HDF5/Zarr backends |
| [`SpatialViz`](@ref) | Makie visualization — pyramid images, scatter, poly, lazy view dispatch |
| `SpatialOmics` | Umbrella — re-exports the full public API |

## Installation

```julia
using Pkg
Pkg.add("SpatialOmics")
```

## Quick start

```julia
using CairoMakie   # load a Makie backend before SpatialViz
using SpatialOmics
import SpatialOmics as SO

# Load from SpatialData OME-ZARR
xen = SO.read(SO.Zarr(), "/path/to/xenium.zarr")

# Inspect structure
keys(xen.images), keys(xen.points), keys(xen.shapes)

# Define a region of interest — lazy, no data copied
roi = view(xen, SpatialExtent(4000, 5000, 1000, 2000))

# Build a composite plot by passing roi to standard Makie verbs
fig = Figure(size=(600, 600))
ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true)
heatmap!(ax, roi.images["morphology_focus"]; channel=1, colormap=:grays)
poly!(ax,    roi.shapes["cell_boundaries"];  color=:transparent, strokecolor=:cyan)
scatter!(ax, roi.points["transcripts"];      markersize=1, color=(:red, 0.3))
fig
```
