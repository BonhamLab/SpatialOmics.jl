# Persist a dataset safely

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
discard!(dataset, "transcripts")

restored = coords(points(dataset, "transcripts"))[1]
close(dataset)
rm(root; recursive=true)
restored
```

For temporary work, `with_dataset() do dataset ... end` guarantees cleanup.
Call `keep!` inside the block when the result should become persistent.

Native stores carry an explicit format version. An unversioned or unsupported
store is rejected with a rebuild instruction; SpatialOmics does not silently
upgrade it. Inspect a store with `native_store_version(path)` before opening it
when its origin is uncertain.
