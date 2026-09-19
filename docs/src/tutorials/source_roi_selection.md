# Select acquisition sources and geometric ROIs

Acquisition identity and geometric containment answer different questions. A
transcript acquired in one FOV does not become an observation from another FOV
merely because their footprints overlap.

## Build two overlapping sources

```@example source-roi-selection
using SpatialOmics

square(xmin, xmax, ymin, ymax) = Polygon(Point2f[
    (xmin, ymin), (xmax, ymin), (xmax, ymax),
    (xmin, ymax), (xmin, ymin),
])

dataset = SpatialDataset()
push!(dataset, CoordinateSystem("global_um"; units=("µm", "µm")))
dataset["fov_footprints"] = SpatialShapes(
    [square(0, 10, 0, 10), square(5, 15, 0, 10)];
    instance_id=Int32[1, 2],
    origins=["fov_1", "fov_2"],
    coord_system="global_um",
)
push!(dataset, AcquisitionSource("fov_1"; region="fov_footprints", instance_id=1))
push!(dataset, AcquisitionSource("fov_2"; region="fov_footprints", instance_id=2))

dataset["transcripts"] = SpatialPoints(
    [Point2f(7, 5), Point2f(7, 5), Point2f(12, 5)];
    origins=["fov_1", "fov_2", "fov_2"],
    coord_system="global_um",
)
nothing
```

The first two transcripts have identical coordinates in the overlap but
different acquisition origins.

## Select exact sources

Index-like source selection returns exactly the requested sources. Selecting
FOVs 1 and 5 would not imply the rectangular region between them.

```@example source-roi-selection
fov_2_transcripts = points(view(dataset, "fov_2"), "transcripts")
(
    coordinates=coords(fov_2_transcripts),
    sources=[source(fov_2_transcripts, i) for i in 1:length(fov_2_transcripts)],
    parent_rows=only(parentindices(fov_2_transcripts)),
)
```

```@example source-roi-selection
both = points(view(dataset, ["fov_1", "fov_2"]), "transcripts")
(length(both), [source(both, i) for i in 1:length(both)])
```

For images and labels, multi-source selection returns
[`SpatialRasterTiles`](@ref): positioned crops with no dense allocation for the
space between distant sources.

## Select by geometry

A geometric ROI ignores acquisition identity and includes every observation
whose current coordinates satisfy the query.

```@example source-roi-selection
overlap = SpatialExtent(6, 8, 4, 6; coord_system="global_um")
overlap_transcripts = points(view(dataset, overlap), "transcripts")
(
    length(overlap_transcripts),
    [source(overlap_transcripts, i) for i in 1:length(overlap_transcripts)],
)
```

Use source selection to answer “what did this acquisition produce?” and an ROI
to answer “what is currently located here?” If an older vector element has no
origin metadata, dataset-level source selection warns and falls back to the
registered footprint; it never silently claims exact provenance.

```@example source-roi-selection
close(dataset; discard=true) # hide
nothing # hide
```
