# Elements

Spatial element types represent the core data modalities in a spatial omics experiment.

## Point collections

```@docs
SpatialPoints
```

### Point accessors

```@docs
coords
features
feature_ids
origins
origin_ids
coord_system
instance_id
instance_ids
with_instance_ids
```

## Shape collections

```@docs
SpatialShapes
SpatialShape
```

### Shape accessors

```@docs
geometries
```

## Utilities

```@docs
subsample
top_features
count_per_instance
```

## Re-exported geometry types

`Polygon`, `MultiPolygon`, and `Point2f` are re-exported from
`GeometryBasics.jl`. Use a `MultiPolygon` when one biological object has
multiple disconnected components; it remains one row in `SpatialShapes`.
