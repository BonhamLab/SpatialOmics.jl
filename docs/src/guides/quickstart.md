# Quickstart

```julia
using SpatialOmics
import SpatialOmics as SO

# Load from SpatialData OME-ZARR
ds = SO.read(SO.Zarr(), "/path/to/experiment.zarr")

# Inspect structure
keys(ds.images), keys(ds.points), keys(ds.tables)

# Export back to SpatialData
SO.write(ds, "/path/to/output.zarr", SO.Zarr())
```
