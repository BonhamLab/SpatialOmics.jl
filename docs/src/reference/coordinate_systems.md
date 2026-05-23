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
translation
scaling
rotation
flip_y
```

## Operations

```@docs
compose
apply
apply!
resolve
```
