# Dataset

The `SpatialDataset` is the root container for a spatial omics experiment.
See [The data model](@ref) for a conceptual overview.

## Container types

```@docs
SpatialDataset
BackingStore
AcquisitionSource
```

## Lifecycle

```@docs
with_dataset
keep!
save!
discard!
edit!
touch!
isdirty
dirty
Base.close(::SpatialDataset)
```

## Accessors

```@docs
elements
coord_systems
transform
sources
source
relations
```

## Typed element retrieval

These functions retrieve a specific named element from a dataset (or dataset
view) and verify its type. Prefer them over `elements(ds)[name]` to catch
element-type mismatches early.

```@docs
images
labels
points
shapes
```
