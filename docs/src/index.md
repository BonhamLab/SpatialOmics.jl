# SpatialOmics.jl

A modular, high-performance Julia ecosystem for spatial transcriptomics analysis.

## Packages

| Package | Role |
|---------|------|
| [`SpatialOmicsBase`](@ref) | Core types — `SpatialDataset`, elements, coordinate systems |
| [`SpatialIO`](@ref) | All I/O: platform readers, SpatialData interop, HDF5/Zarr backends |
| `SpatialOmics` | Umbrella — re-exports the full public API |

## Installation

```julia
using Pkg
Pkg.add("SpatialOmics")
```

## Quick start

```julia
using SpatialOmics

# Load a CosMx experiment
ds = load_cosmx("/path/to/DecodedFiles_run/")

# Load from SpatialData OME-ZARR
ds = from_spatialdata("/path/to/experiment.zarr")

# Export to SpatialData OME-ZARR
to_spatialdata(ds, "/path/to/output.zarr")
```
