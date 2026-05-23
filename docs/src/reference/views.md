# Views and extents

See [Lazy views and spatial filtering](@ref) for a conceptual overview.

## Extents and regions

```@docs
SpatialExtent
SpatialROI
```

## Views

```@docs
SpatialDatasetView
SpatialElementView
```

## Accessors

```@docs
geometry
select
```

## Example: constructing and intersecting extents

```@example views
using SpatialOmics

a = SpatialExtent(0.0, 1000.0, 0.0, 500.0; coord_system="global")
b = SpatialExtent(500.0, 1500.0, 0.0, 500.0; coord_system="global")

intersect(a, b)   # SpatialExtent(500.0, 1000.0, 0.0, 500.0, "global")
```

```@example views
union(a, b)       # SpatialExtent(0.0, 1500.0, 0.0, 500.0, "global")
```
