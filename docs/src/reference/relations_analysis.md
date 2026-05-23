# Relations and analysis

`SpatialRelation` is the generic structure for weighted associations between
spatial elements. The `RelationKind` dispatch tokens select the algorithm used
by `analyze`.

## Relation kinds

```@docs
RelationKind
Membership
Proximity
KNN
Expression
```

## Relations

```@docs
SpatialRelation
nobs
nvar
var_names
annotate
```

## Analysis

```@docs
analyze
distances
PointDensity
density
ShapeColorView
```
