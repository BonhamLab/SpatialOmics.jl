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
| `image!(ax, channel(img, 1))` | Pyramid-aware image; selects the correct level on each zoom |
| `scatter!(ax, pts)` | Transcript / centroid scatter from `coords(pts)` |
| `poly!(ax, shp)` | Cell boundary / annotation polygons from `SpatialShapes` |
| `heatmap!(ax, lbl)` | Integer label overlay from `SpatialLabels` |

All verbs have mutating `!` variants and accept the full set of Makie keyword
arguments (`color`, `strokewidth`, `markersize`, etc.).

## Pyramid-aware images

`SpatialImage` objects loaded from OME-Zarr carry pre-computed pyramid levels
as lazy `DiskArray`-backed arrays. `image!(ax, img)` selects the correct
resolution level on every zoom or pan event. No pixels are loaded until
a viewport is established.

```julia
img = images(ds, "morphology_focus")

fig, ax, _ = image(channel(img, 1);
    axis=(; aspect=DataAspect(), yreversed=true),
    figure=(; size=(700, 230)))
tightlimits!(ax)
```

## Channel selection and display scaling

For multi-channel images, extract a single channel with `channel`, then apply
`scaleminmax` for display-time intensity normalisation:

```julia
dapi = channel(img, 1)           # or channel(img, "DAPI")
image!(ax, scaleminmax(dapi))
```

`scaleminmax` samples the intensity range from the coarsest pyramid level and
attaches a min-max display transform applied at render time — no copy is made.

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
image!(ax,   scaleminmax(channel(images(roi, "morphology_focus"), 1)))
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

fig = Figure(size=(900, 250))
for (i, lbl) in enumerate(ch_names)
    ax = Axis(fig[1, i]; title=lbl, aspect=DataAspect(),
              yreversed=true, xticksvisible=false, yticksvisible=false)
    image!(ax, scaleminmax(channel(images(ds, "morphology_focus"), i)))
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
image!(ax, colorview(RGB, r, g, b))
```
