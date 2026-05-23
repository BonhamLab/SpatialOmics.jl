# The data model

## Why a dataset container?

Spatial omics data is not a single array. A typical Xenium or CosMx experiment
produces transcripts (millions of 2-D points), cell boundaries (tens of thousands
of polygons), a morphology image (gigapixels), a segmentation mask, and a cell × gene
count matrix — each at different resolutions, in different coordinate spaces, with
cross-referencing instance IDs. Any representation that flattens this into one
structure either loses information or forces you to carry alignment metadata by hand.

`SpatialDataset` is an explicit container that holds all of these together, named
and typed, without prejudicing any one modality as primary.

## Element types

| Type | What it represents |
|------|--------------------|
| `SpatialPoints` | Point cloud — transcripts, cell centroids, any 2-D coordinates with optional feature labels and instance IDs |
| `SpatialShapes` | Polygon collection — cell boundaries, tissue annotations, FOV outlines |
| `SpatialImage` | Raster image — morphology, DAPI, IF channels; optionally multi-channel and pyramid-backed |
| `SpatialLabels` | Integer segmentation mask — pixel-level cell assignment, same coordinate space as a paired image |

All four types belong to one named coordinate system. A dataset can hold any
number of elements of each type, keyed by name.

## Typed accessors

Elements are stored internally in a flat `OrderedDict`. The typed accessor
functions — `points(ds, name)`, `shapes(ds, name)`, `images(ds, name)`,
`labels(ds, name)` — retrieve a specific element and verify its type, raising
a clear error if the name maps to the wrong kind:

```julia
tx  = points(ds, "transcripts")    # SpatialPoints or error
img = images(ds, "morphology")     # SpatialImage or error
```

Accessing `elements(ds)` directly returns the raw `OrderedDict` without type
checking. Use typed accessors in application code; `elements` is useful for
iteration or introspection.

## The backing store

Every `SpatialDataset` is associated with a `BackingStore` — a Zarr directory
on disk. This is not optional. The design exists because:

1. **Lazy loading** — images and large point clouds can exceed available RAM.
   Zarr arrays are read on demand through `DiskArrays.jl`.
2. **Persistence by default** — operations that produce new datasets (such as
   `read`) always have a place to write without a separate "save" step.
3. **SpatialData compatibility** — the on-disk layout matches the
   [SpatialData specification](https://spatialdata.scverse.org/), enabling
   round-trip with Python tools without a conversion step.

When `SpatialDataset()` is called without a `path`, a temporary directory is
created and owned by the dataset — it is deleted automatically when the dataset
is garbage collected or `close`d. Supply `path` to write directly to a
persistent location, or call `write!(ds, path, SpatialDataZarr())` to move a
temporary store to a permanent one.

## Instance IDs and cross-element linkage

The `instance_id` field threads through all element types. A cell segmented in
`SpatialLabels` has a pixel label; `SpatialShapes` carries the same value as
`instance_id` for the cell boundary polygon; `SpatialPoints` assigns transcripts
to cells via `instance_id`. `SpatialRelation` aggregates these into a count
matrix keyed by the same IDs.

This design avoids duplicate data: shapes, labels, transcripts, and count
matrices all refer to the same cells by a shared integer key rather than
embedding coordinates or names redundantly.
