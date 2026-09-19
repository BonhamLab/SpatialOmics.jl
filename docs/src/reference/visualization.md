# Visualization types

These types and functions support image display. Makie plot verbs (`heatmap!`,
`scatter!`, `poly!`) are defined in the Makie extension and dispatch on
`SpatialImage`, `SpatialPoints`, and `SpatialShapes` directly — see the
[Visualization guide](../guides/visualization.md) for usage patterns.

## Image types

```@docs
SpatialImage
SpatialLabels
SpatialRasterTiles
SpatialImageColorView
```

## Image utilities

```@docs
nchannels
channel_names
channel
scaleminmax
colorview
build_pyramid!
ensure_pyramid!
data
```

## Re-exported color types

`Gray` and `RGB` are re-exported from `Colors.jl`. Use them as the colorant
argument to `colorview`:

```julia
cview = colorview(Gray, scaleminmax(channel(img, 1)))
rgb   = colorview(RGB,  channel(img, 1), channel(img, 2), channel(img, 3))
```
