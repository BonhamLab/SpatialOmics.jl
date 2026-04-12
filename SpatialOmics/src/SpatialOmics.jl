"""
    SpatialOmics

User-facing umbrella package for the spatial transcriptomics ecosystem.

## Two-tier API

**Tier 1 — types and data manipulation** (`using SpatialOmics`):
Core types and dataset mutation functions are exported and land directly in
scope.

**Tier 2 — I/O** (`import SpatialOmics as SO`):
All I/O is available as qualified names only — nothing in tier 2 is exported.

```julia
using SpatialOmics          # SpatialDataset, images, points, ds["key"] = el, …
import SpatialOmics as SO   # SO.load_cosmx, SO.read, SO.write, SO.Zarr, …

ds = SO.load_cosmx("/data/run/")
SO.write(ds, "/out/experiment.zarr", SO.Zarr())
ds2 = SO.read(SO.Zarr(), "/out/experiment.zarr")
```
"""
module SpatialOmics

# ── Core data structures (exported) ──────────────────────────────────────────

import SpatialOmicsBase:
    SpatialDataset,
    SpatialElement,
    SpatialImage, SpatialPoints, SpatialLabels, SpatialShapes, SpatialTable

export SpatialDataset,
       SpatialElement,
       SpatialImage, SpatialPoints, SpatialLabels, SpatialShapes, SpatialTable

# ── Dataset construction and mutation (exported) ──────────────────────────────

import SpatialOmicsBase:
    spatial_dataset,
    images, labels, points, shapes, tables, metadata,
    channels, channels!

export spatial_dataset,
       images, labels, points, shapes, tables, metadata,
       channels, channels!

# ── Coordinate systems and transformations (exported) ─────────────────────────

import SpatialOmicsBase:
    CoordinateSystem, Transformation,
    AffineTransformation, IdentityTransformation,
    transform_coordinates, compose_transformations, invert_transformation

export CoordinateSystem, Transformation,
       AffineTransformation, IdentityTransformation,
       transform_coordinates, compose_transformations, invert_transformation

# ── Spatial query primitives (exported) ───────────────────────────────────────

import SpatialOmicsBase:
    SpatialExtent,
    SpatialElementView, SpatialDatasetView,
    extent, intersects,
    region, region_key, instance_key,
    ImagePyramidSampler

export SpatialExtent,
       SpatialElementView, SpatialDatasetView,
       extent, intersects,
       region, region_key, instance_key,
       ImagePyramidSampler

# ── Tier 2 — spatial operations (unexported/public) ──────────────────────────
# Use qualified: SO.filter(roi, key)
import SpatialOmicsBase: filter
public filter

# ── Tier 2 — I/O platform loaders (unexported/public) ────────────────────────
import SpatialIO:
    PlatformReader,
    XeniumReader, VisiumReader, CosMxReader, MerfishReader,
    validate_path, read_data, load_spatial_data,
    load_xenium, load_visium, load_cosmx, load_merfish

public PlatformReader
public XeniumReader, VisiumReader, CosMxReader, MerfishReader
public validate_path, read_data, load_spatial_data
public load_xenium, load_visium, load_cosmx, load_merfish

# ── Tier 2 — I/O format tokens (unexported/public) ───────────────────────────
import SpatialIO:
    StorageFormat,
    Zarr, NativeH5, AnnData, GeoJSON,
    TranscriptsParquet, CellsParquet, CellFeatureMatrixH5, TissuePositionsCSV

public StorageFormat
public Zarr, NativeH5, AnnData, GeoJSON
public TranscriptsParquet, CellsParquet, CellFeatureMatrixH5, TissuePositionsCSV

# ── Tier 2 — I/O read/write (unexported/public) ──────────────────────────────
import SpatialIO: read, write
public read, write

end # module SpatialOmics
