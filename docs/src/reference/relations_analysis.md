# Relations and analysis

`SpatialRelation` is the generic structure for weighted associations between
spatial elements. The `RelationKind` dispatch tokens select the algorithm used
by `analyze`.

## Relation kinds

```@docs
RelationKind
Membership
Expression
```

## Relations

```@docs
SpatialRelation
source_ids
destination_ids
nobs
nvar
obs_names
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
