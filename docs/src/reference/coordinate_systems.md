# Coordinate systems

See [Coordinate systems](../explanation/coordinate_systems.md) for a conceptual overview of the transform graph design.

## Coordinate system

```@docs
CoordinateSystem
```

## Transformation types

```@docs
AbstractTransformation
Identity
Affine
Sequence
```

## Transformation constructors

```@docs
SpatialOmics.translation
SpatialOmics.scaling
SpatialOmics.rotation
SpatialOmics.flip_y
```

## Operations

```@docs
SpatialOmics.compose
apply
apply!
resolve
```
