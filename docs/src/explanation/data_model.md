# The data model

## Why a dataset container?

Spatial omics data is not a single array. A typical Xenium or CosMx experiment
produces transcripts (millions of 2-D points), cell boundaries (tens of thousands
of polygons), a tissue image (gigapixels), a segmentation mask, and a cell × gene
count matrix — each at different resolutions, in different coordinate spaces, with
cross-referencing instance IDs. Any representation that flattens this into one
structure either loses information or forces you to carry alignment metadata by hand.

[`SpatialDataset`](@ref) is an explicit container that holds all of these together, named
and typed, without prejudicing any one modality as primary.

## Element types

| Type | What it represents |
|------|--------------------|
| [`SpatialPoints`](@ref) | Point cloud — transcripts, cell centroids, any 2-D coordinates with optional feature labels and instance IDs |
| [`SpatialShapes`](@ref) | Polygon collection — cell boundaries, tissue annotations, FOV outlines |
| [`SpatialImage`](@ref) | Raster image — tissue images, fluorescence channels; optionally multi-channel and pyramid-backed |
| [`SpatialLabels`](@ref) | Integer segmentation mask — pixel-level cell assignment, same coordinate space as a paired image |

All four types belong to one named coordinate system. A dataset can hold any
number of elements of each type, keyed by name.

## Typed accessors

Elements are stored internally in a flat `OrderedDict`. The typed accessor
functions — [`points`](@ref), [`shapes`](@ref), [`images`](@ref),
[`labels`](@ref) — retrieve a specific element and verify its type, raising
a clear error if the name maps to the wrong kind:

```julia
tx  = points(ds, "transcripts")    # SpatialPoints or error
img = images(ds, "morphology")     # SpatialImage or error
```

[`elements`](@ref)`(ds)` returns a shallow dictionary snapshot. Changing that
dictionary does not change the dataset. Use typed accessors for retrieval and
[`edit!`](@ref) when mutating an attached element.

## The backing store

Every [`SpatialDataset`](@ref) is associated with a [`BackingStore`](@ref) — a
Zarr directory on disk. This provides a durable target for large, lazily loaded
data without forcing every analysis step to perform I/O.

1. **Lazy loading** — images and large point clouds can exceed available RAM.
   Zarr arrays are read on demand through `DiskArrays.jl`.
2. **Explicit checkpoints** — supported mutations are staged and can be saved
   together or element by element.
3. **Visible state** — [`isdirty`](@ref), [`dirty`](@ref), and dataset display
   distinguish saved data from unsaved work.

When `SpatialDataset()` is called without a `path`, a temporary directory is
created and owned by the dataset. Supply `path` to choose a persistent backing
location. In both cases, mutations remain staged until [`save!`](@ref) is
called:

```julia
ds["transcripts"] = transcripts
isdirty(ds)                  # true
dirty(ds)                    # identifies the staged element
save!(ds, "transcripts")     # save one element
save!(ds)                    # save everything else
```

Use `save!(ds; path="/data/experiment.zarr")` to atomically write a complete
snapshot and rebind a temporary dataset to a permanent location. `close(ds)`
rejects unsaved changes; `discard!(ds)` restores saved state, while
`close(ds; discard=true)` explicitly abandons it.

Mutation through package operations is tracked automatically. For mutation
through an external API, use a scoped edit or mark the element afterward:

```julia
edit!(ds, "transcripts") do points
    points.feature_id[1] = 2
end

external_mutation!(points(ds, "transcripts"))
touch!(ds, "transcripts")
```

The native layout is a SpatialOmics format. Reading supported Python
SpatialData stores is an import operation; native stores should not be assumed
to round-trip through Python without an explicit exporter.

## Instance IDs and cross-element linkage

The [`instance_id`](@ref) field threads through all element types. A cell segmented in
[`SpatialLabels`](@ref) has a pixel label; [`SpatialShapes`](@ref) carries the same value as
[`instance_id`](@ref) for the cell boundary polygon; [`SpatialPoints`](@ref) assigns transcripts
to cells via [`instance_id`](@ref). [`SpatialRelation`](@ref) aggregates these into a count
matrix keyed by the same IDs.

This design avoids duplicate data: shapes, labels, transcripts, and count
matrices all refer to the same cells by a shared integer key rather than
embedding coordinates or names redundantly.
