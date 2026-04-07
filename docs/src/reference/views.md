# Spatial views and extents

## Bounding boxes

```@docs
SpatialExtent
intersects
extent
```

## Lazy views

`view(ds, ext)` and `view(el, ext)` return lazy views analogous to Julia's
`SubArray` from `@view x[...]`. No data is copied; the spatial filter is
applied only when a plot verb or `crop` materialises the underlying array.

```@docs
SpatialDatasetView
SpatialElementView
```

## Materialisation

`crop` returns a concrete filtered copy. Use it when you need a standalone
element outside of a plotting context.

```@docs
crop
```
