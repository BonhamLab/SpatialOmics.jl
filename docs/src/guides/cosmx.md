# CosMx workflow

```julia
using SpatialOmics
import SpatialOmics as SO

# Load a CosMx DecodedFiles export
ds = SO.load_cosmx("/path/to/run_dir/")

# Explore transcripts for a single FOV
tx  = ds.points["transcripts"]
fov = tx.features[tx.features.fov .== 1, :]

# Export to SpatialData OME-ZARR
SO.write(ds, "/path/to/output.zarr", SO.Zarr())
```
