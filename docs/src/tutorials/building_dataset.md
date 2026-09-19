# Build a spatial dataset

This tutorial builds a small transcript-and-cell dataset using only public
constructors and accessors. The same pattern is useful for an unsupported assay
format: parse the source tables at the boundary, then construct ordinary
SpatialOmics elements.

## Construct transcript points

`SpatialPoints` accepts any Tables.jl-compatible object. Identify the feature
column with `gene` and keep other per-transcript columns in `features`.
Acquisition origins are independent of spatial coordinates.

```@example building-dataset
using SpatialOmics

transcript_table = (
    x = Float32[1, 2, 7, 8],
    y = Float32[1, 2, 1, 2],
    gene = ["Actb", "Gapdh", "Actb", "Krt8"],
    quality = Float32[0.98, 0.93, 0.96, 0.91],
    origin = fill("fov_1", 4),
)

transcripts = SpatialPoints(
    transcript_table;
    gene=:gene,
    origin=:origin,
    features=(quality=transcript_table.quality,),
    coord_system="global_um",
)

(length(transcripts), features(transcripts), features(transcripts, :quality))
```

Feature labels use a compact codebook plus integer IDs. Use `features` and
`feature_ids` instead of depending on those storage fields directly.

## Construct cell polygons

`Polygon` and `MultiPolygon` are re-exported from GeometryBasics. One entry in
`SpatialShapes` represents one object, even when that object is a
`MultiPolygon` with disconnected components.

```@example building-dataset
square(xmin, xmax, ymin, ymax) = Polygon(Point2f[
    (xmin, ymin), (xmax, ymin), (xmax, ymax),
    (xmin, ymax), (xmin, ymin),
])

cells = SpatialShapes(
    [square(0, 4, 0, 4), square(6, 10, 0, 4)];
    instance_id=Int32[101, 102],
    origins=["fov_1", "fov_1"],
    coord_system="global_um",
)

(length(cells), instance_id(cells), coord_system(cells))
```

## Assemble the dataset

Register coordinate systems and acquisition sources explicitly. The source
footprint is a spatial object; the source name is the acquisition identity
stored on observations.

```@example building-dataset
dataset = SpatialDataset(metadata=Dict("technology" => "synthetic example"))
push!(dataset, CoordinateSystem("global_um"; units=("µm", "µm")))

dataset["transcripts"] = transcripts
dataset["cells"] = cells
dataset["fov_footprints"] = SpatialShapes(
    [square(0, 10, 0, 4)];
    instance_id=Int32[1],
    origins=["fov_1"],
    coord_system="global_um",
)
push!(dataset, AcquisitionSource(
    "fov_1";
    region="fov_footprints",
    instance_id=1,
    attributes=Dict("vendor_fov" => 1),
))

(
    element_names=collect(keys(elements(dataset))),
    coordinate_systems=coord_systems(dataset),
    acquisition_sources=sources(dataset),
)
```

Typed accessors fail early if a name refers to the wrong kind of element:

```@example building-dataset
tx = points(dataset, "transcripts")
bounds = shapes(dataset, "cells")
(top_features(tx, 3), count_per_instance(tx))
```

These transcript rows have not yet been assigned to cells, so their
`instance_id` values are zero and `count_per_instance` is empty. Spatial
assignment is covered in [Assign transcripts and summarise expression](@ref).

```@example building-dataset
close(dataset; discard=true) # hide
nothing # hide
```
