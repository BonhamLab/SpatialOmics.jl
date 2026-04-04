# SpatialData / OME-ZARR I/O

Format-level persistence uses `SpatialOmics.read` / `SpatialOmics.write` dispatched
on format type tokens.  These functions are **not** exported — use
module-qualified access:

```julia
import SpatialOmics as SO

# Read from OME-ZARR (SpatialData format)
ds = SO.read(SO.Zarr(), "/data/experiment.zarr")

# Read from native SpatialIO HDF5
ds = SO.read(SO.NativeH5(), "experiment.h5")

# Write to OME-ZARR
SO.write(ds, "/out/experiment.zarr", SO.Zarr())

# Write to native HDF5
SO.write(ds, "experiment.h5", SO.NativeH5())

# Export to AnnData .h5ad
SO.write(ds, "table.h5ad", SO.AnnData())
```

## Format type tokens

```@docs
SpatialOmics.Zarr
SpatialOmics.NativeH5
```

## Read / write

```@docs
SpatialOmics.read
SpatialOmics.write
```
