# Persist a dataset safely

```@setup persistence
using CairoMakie
using Markdown
CairoMakie.activate!(type="svg")
set_theme!(Theme(
    fontsize=15,
    Figure=(; backgroundcolor=:white),
    Axis=(; xgridvisible=false, ygridvisible=false),
))
```

Every `SpatialDataset` has a backing store, but in-memory changes are staged
until `save!` is called. Closing a dirty dataset refuses to discard those
changes implicitly.

## Save and reopen

Use a persistent path when constructing the dataset or provide one to `save!`.

```@example persistence
using SpatialOmics

root = mktempdir()
path = joinpath(root, "example.zarr")
dataset = SpatialDataset(path=path)
push!(dataset, CoordinateSystem("global_um"; units=("µm", "µm")))
dataset["transcripts"] = SpatialPoints(
    [Point2f(1, 2), Point2f(3, 4)];
    coord_system="global_um",
)
dataset.metadata["sample"] = "example"

before_save = isdirty(dataset)
save!(dataset)
after_save = isdirty(dataset)
close(dataset)

reopened = read(SpatialDataZarr(), path)
result = (
    before_save=before_save,
    after_save=after_save,
    sample=reopened.metadata["sample"],
    transcripts=length(points(reopened, "transcripts")),
)
close(reopened)
result
```

## Track edits explicitly

Package-mediated operations mark attached elements dirty. For a direct edit,
use `edit!` so the dataset knows that the element must be saved.

```@example persistence
dataset = read(SpatialDataZarr(), path)
edit!(dataset, "transcripts") do transcripts
    coords(transcripts)[1] = Point2f(9, 9)
end

changes = [(change.kind, change.name, change.state) for change in dirty(dataset)]
save!(dataset, "transcripts")
saved_coordinates = copy(coords(points(dataset, "transcripts")))
(changes=changes, dirty_after_save=isdirty(dataset))
```

`save!(dataset, "transcripts")` persists one named element. Other staged
changes, if any, remain dirty. `save!(dataset)` saves the complete pending
change set.

## Discard deliberately

`discard!` reloads saved state. `close(dataset; discard=true)` abandons all
remaining changes and is most useful for temporary exploratory datasets.

```@example persistence
edit!(dataset, "transcripts") do transcripts
    coords(transcripts)[1] = Point2f(99, 99)
end
staged_coordinates = copy(coords(points(dataset, "transcripts")))
discard!(dataset, "transcripts")

restored_coordinates = copy(coords(points(dataset, "transcripts")))
restored = restored_coordinates[1]
close(dataset)
rm(root; recursive=true)
restored
```

`discard!` restores the saved element rather than silently retaining the staged
edit. All three panels use the same limits so the discarded displacement is
apparent.

```@eval persistence
figure = Figure(size=(900, 310))
states = (("saved", saved_coordinates, :seagreen),
          ("staged edit", staged_coordinates, :darkorange),
          ("after discard!", restored_coordinates, :steelblue))
for (column, (title, coordinates, color)) in enumerate(states)
    axis = Axis(figure[1, column]; aspect=DataAspect(), title,
                xlabel="x", ylabel=column == 1 ? "y" : "")
    scatter!(axis, coordinates[2:end]; color=:gray60, markersize=12)
    scatter!(axis, coordinates[1:1]; color, markersize=16)
    point = first(coordinates)
    near_edge = point[1] > 90 || point[2] > 90
    offset = near_edge ? -3 : 3
    alignment = near_edge ? (:right, :top) : (:left, :bottom)
    text!(axis, point[1] + offset, point[2] + offset;
          text="($(Int(point[1])), $(Int(point[2])))",
          align=alignment, fontsize=12)
    xlims!(axis, 0, 105); ylims!(axis, 0, 105)
end
save("persistence-states.svg", figure)
Markdown.parse("![Saved, staged, and restored coordinates](persistence-states.svg)")
```

For temporary work, `with_dataset() do dataset ... end` guarantees cleanup.
Call `keep!` inside the block when the result should become persistent.

Native stores carry an explicit format version. An unversioned or unsupported
store is rejected with a rebuild instruction; SpatialOmics does not silently
upgrade it. Inspect a store with `native_store_version(path)` before opening it
when its origin is uncertain.
