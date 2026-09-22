# Select acquisition sources and geometric ROIs

```@setup source-roi-selection
using CairoMakie
using Markdown
CairoMakie.activate!(type="svg")
set_theme!(Theme(
    fontsize=15,
    Figure=(; backgroundcolor=:white),
    Axis=(; xgridvisible=false, ygridvisible=false),
))
```

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

The plot makes the distinction visible: source membership is attached to each
observation, while the ROI is a geometric query over the overlap.

```@eval source-roi-selection
figure = Figure(size=(760, 390))
axis = Axis(figure[1, 1]; aspect=DataAspect(), xlabel="x (µm)", ylabel="y (µm)")
footprints = geometries(shapes(dataset, "fov_footprints"))
poly!(axis, [footprints[1]]; color=(:dodgerblue, 0.18),
      strokecolor=:dodgerblue3, strokewidth=2, label="FOV 1 footprint")
poly!(axis, [footprints[2]]; color=(:darkorange, 0.18),
      strokecolor=:darkorange3, strokewidth=2, label="FOV 2 footprint")
scatter!(axis, [Point2f(7, 5)]; color=:dodgerblue3, marker=:circle,
         markersize=18, label="transcript from FOV 1")
scatter!(axis, [Point2f(7, 5), Point2f(12, 5)]; color=:darkorange3,
         marker=:xcross, markersize=20, label="transcript from FOV 2")
poly!(axis, [Rect2f(6, 4, 2, 2)]; color=(:purple, 0.08),
      strokecolor=:purple, strokewidth=3, linestyle=:dash, label="geometric ROI")
axislegend(axis; position=:rt, framevisible=false, labelsize=12)
xlims!(axis, -0.5, 15.5); ylims!(axis, -0.5, 10.5)
save("source-roi-selection.svg", figure)
Markdown.parse("![Acquisition footprints and geometric ROI](source-roi-selection.svg)")
```

Use source selection to answer “what did this acquisition produce?” and an ROI
to answer “what is currently located here?” If an older vector element has no
origin metadata, dataset-level source selection warns and falls back to the
registered footprint; it never silently claims exact provenance.

```@example source-roi-selection
close(dataset; discard=true) # hide
nothing # hide
```
