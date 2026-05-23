# Visualization

`SpatialOmics` extends standard Makie verbs to dispatch directly on spatial
element types. Users write the same Makie code they already know; the type
system routes to spatial-aware implementations.

## Backend choice

`SpatialOmics` depends only on the abstract `Makie` package. Load a concrete
backend **before** calling any plot verb:

```julia
using GLMakie      # recommended for interactive exploration
# or
using CairoMakie   # recommended for notebooks / publication figures
# or
using WGLMakie     # for Pluto / Jupyter
```

## Plotting verbs

| Call | What it does |
|------|-------------|
| `heatmap!(ax, img; channel=1)` | Pyramid-aware image; `Makie.Resampler` selects the correct level on each zoom |
| `scatter!(ax, pts)` | Transcript / centroid scatter from `coords(pts)` |
| `poly!(ax, shp)` | Cell boundary / annotation polygons from `SpatialShapes` |
| `heatmap!(ax, lbl)` | Integer label overlay from `SpatialLabels` |

All verbs have mutating `!` variants and accept the full set of Makie keyword
arguments (`colormap`, `color`, `strokewidth`, `markersize`, etc.).

## Pyramid-aware images

`SpatialImage` objects loaded from OME-Zarr carry pre-computed pyramid levels
as lazy `DiskArray`-backed arrays. `heatmap!(ax, img)` wraps them in an
`ImagePyramidSampler` and passes it to `Makie.Resampler`, which selects the
correct resolution level on every zoom or pan event. No pixels are loaded until
a viewport is established.

```julia
img = images(ds, "morphology_focus")

fig, ax, _ = heatmap(img; channel=1, colormap=:grays,
    axis=(; aspect=DataAspect(), yreversed=true),
    figure=(; size=(700, 230)))
tightlimits!(ax)
```

## Channel selection

For multi-channel images, pass `channel=i` (integer index) or
`channel="DAPI"` (channel name) to the plot verb. Alternatively, use the
`channel` function to extract a 2-D single-channel `SpatialImage` first:

```julia
dapi = channel(img, 1)           # or channel(img, "DAPI")
heatmap!(ax, scaleminmax(dapi); colormap=:grays)
```

`scaleminmax` samples the intensity range from the coarsest pyramid level and
attaches a min-max display transform that is applied at render time.

## Lazy spatial views

`view(ds, ext)` returns a `SpatialDatasetView` scoped to a bounding box. This
is the preferred way to build composite panels: define the ROI once, then pass
the view to each layer independently.

```julia
ext = SpatialExtent(cx - 500.0, cx + 500.0, cy - 500.0, cy + 500.0; coord_system="global")
roi = view(ds, ext)

fig = Figure(size=(600, 600))
ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true,
           xlabel="x (µm)", ylabel="y (µm)")

# Each call applies the crop independently — no intermediate copies
heatmap!(ax, images(roi, "morphology_focus"); channel=1, colormap=:grays)
poly!(ax,    shapes(roi, "cell_boundaries");  color=:transparent,
                                              strokecolor=:cyan, strokewidth=0.4)
scatter!(ax, points(roi, "transcripts");      color=(:red, 0.25), markersize=1)
tightlimits!(ax)
fig
```

When a `SpatialImage` is accessed through a `SpatialDatasetView`, the axis
limits are automatically constrained to the view extent.

## Multi-channel display

```julia
ch_names = ["DAPI", "ATP1A1/CD45/E-Cad", "18S", "AlphaSMA/Vim"]
cmaps    = [:grays, :viridis, :magma, :plasma]

fig = Figure(size=(900, 250))
for (i, (lbl, cmap)) in enumerate(zip(ch_names, cmaps))
    ax = Axis(fig[1, i]; title=lbl, aspect=DataAspect(),
              yreversed=true, xticksvisible=false, yticksvisible=false)
    heatmap!(ax, images(ds, "morphology_focus"); channel=i, colormap=cmap)
    tightlimits!(ax)
end
fig
```

## RGB composite

Use `colorview` to combine single-channel images into a colour composite:

```julia
r = scaleminmax(channel(img, 1))
g = scaleminmax(channel(img, 2))
b = scaleminmax(channel(img, 3))
heatmap!(ax, colorview(RGB, r, g, b))
```
