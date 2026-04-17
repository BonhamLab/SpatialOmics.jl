# CosMx workflow

```julia
using SpatialOmics
import SpatialOmics as SO

# Load a CosMx raw export
ds = SO.load(SO.CosMxReader(), "/path/to/cosmx_export/")

# Include Morphology2D TIF tiles (written to zarr on SO.write)
ds = SO.load(SO.CosMxReader(; morphology_dir="/path/to/Morphology2D"),
             "/path/to/cosmx_export/")

# Explore transcripts
tx = ds.points["transcripts"]

# Export to SpatialData OME-ZARR
SO.write(ds, "/path/to/output.zarr", SO.Zarr())
```

For datasets without a dedicated reader (e.g. custom segmentation output, or
the CosMx per-FOV DecodedFiles hierarchy), build a `SpatialDataset` directly
from CSV files — see the [manual ingest walkthrough](@ref).
