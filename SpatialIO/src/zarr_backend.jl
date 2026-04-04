# SpatialIO/src/zarr_backend.jl
#
# open_zarr / write_zarr — thin wrappers around the SpatialData native backend.
#
# The on-disk format is the SpatialData OME-ZARR layout (Zarr v3 with SpatialData
# metadata conventions). These functions exist so Julia code can use the intuitive
# open_zarr / write_zarr names without having to know the SpatialData origin.

"""
    open_zarr(path::String; use_python::Bool=false) -> SpatialDataset

Open a SpatialData-compatible OME-ZARR store at `path` and return a
`SpatialDataset`.  Large arrays (images, labels) are backed by lazy
`ZarrV3Array` objects and are not loaded into memory until accessed.

When `use_python=true` the `spatialdata` Python package is used via PythonCall
for full spec compliance.  The default native backend covers all elements
written by `write_zarr` / `to_spatialdata`.

# Example
```julia
ds = open_zarr("/data/my_experiment.zarr")
```
"""
open_zarr(path::String; kwargs...) = from_spatialdata(path; kwargs...)

"""
    write_zarr(ds::SpatialDataset, path::String; use_python::Bool=false)

Serialise `ds` to a SpatialData-compatible OME-ZARR store at `path`.

# Example
```julia
write_zarr(ds, "/output/my_experiment.zarr")
```
"""
write_zarr(ds::SpatialDataset, path::String; kwargs...) =
    to_spatialdata(ds, path; kwargs...)
