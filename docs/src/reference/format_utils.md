# Format utilities

Format-level I/O uses `SpatialOmics.read` / `SpatialOmics.write` dispatched on
format type tokens.  These functions are **not** exported — use module-qualified
access:

```julia
import SpatialOmics as SO

# Read a transcripts parquet file
tx = SO.read(SO.TranscriptsParquet(), "/path/to/transcripts.parquet")

# Read a Visium tissue positions CSV
pos = SO.read(SO.TissuePositionsCSV(), "/path/to/tissue_positions.csv")

# Read a 10x cell-feature-barcode HDF5
matrix, barcodes, features = SO.read(SO.CellFeatureMatrixH5(), "matrix.h5")
```

## Format type tokens

```@docs
SpatialOmics.StorageFormat
SpatialOmics.TranscriptsParquet
SpatialOmics.CellsParquet
SpatialOmics.CellFeatureMatrixH5
SpatialOmics.TissuePositionsCSV
SpatialOmics.AnnData
SpatialOmics.GeoJSON
```
