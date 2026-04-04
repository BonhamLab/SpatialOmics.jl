# Data structures

## Dataset

```@docs
SpatialDataset
spatial_dataset
add_image!
add_labels!
add_points!
add_shapes!
add_table!
get_image
get_labels
get_points
get_shapes
get_table
get_metadata
set_metadata!
```

## Elements

```@docs
SpatialElement
SpatialImage
SpatialPoints
SpatialLabels
SpatialShapes
SpatialTable
```

## Coordinate systems and transformations

```@docs
CoordinateSystem
Transformation
AffineTransformation
IdentityTransformation
transform_coordinates
compose_transformations
invert_transformation
```
