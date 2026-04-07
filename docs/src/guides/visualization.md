# Visualization

`SpatialViz.jl` extends standard Makie verbs to dispatch directly on
SpatialOmicsBase element types. Users write the same Makie code they already
know; the type system routes to spatial-aware implementations.

## Backend choice

`SpatialViz` depends only on the abstract `Makie` package. Load a concrete
backend **before** `using SpatialViz`:

```julia
using GLMakie      # recommended for interactive exploration
# or
using CairoMakie   # recommended for notebooks / publication figures
# or
using WGLMakie     # for Pluto / Jupyter

using SpatialViz
```

## Plotting verbs

| Call | What it does |
|------|-------------|
| `heatmap(img; channel=1)` | Pyramid-aware image; `Makie.Resampler` selects the correct level on each zoom |
| `scatter(pts)` | Transcript / centroid scatter from first two coordinate columns |
| `poly(shp)` | Cell boundary / annotation polygons from `SpatialShapes` |
| `heatmap(lbl)` | Integer label overlay from `SpatialLabels` |

All verbs have mutating `!` variants and accept the full set of Makie keyword
arguments (`colormap`, `color`, `strokewidth`, `markersize`, etc.).

## Pyramid-aware images

`SpatialImage` objects loaded from OME-Zarr carry pre-computed pyramid levels
as lazy `DiskArray`-backed arrays. `heatmap(img)` wraps them in
`ImagePyramidSampler` and passes the sampler to `Makie.Resampler`, which
selects the correct resolution level on every zoom or pan event. No pixels
are loaded until a viewport is established.

```julia
img = xen.images["morphology_focus"]

# Channel 1 (default)
fig, ax, plt = heatmap(img; colormap=:grays,
    axis=(; aspect=DataAspect()), figure=(; size=(700, 230)))

# Select a different channel
heatmap!(ax, img; channel=2, colormap=:viridis)
```

## Lazy spatial views

`view(ds, extent)` returns a `SpatialDatasetView` — a lazy pointer scoped to
a bounding box. Accessing `.images["key"]`, `.points["key"]`, etc. returns a
`SpatialElementView{T}` that carries the parent element and the extent. The
spatial filter is applied only when a plot verb materialises the data.

This is the preferred way to build composite panels: define the ROI once, then
pass it to each layer independently.

```julia
roi = view(xen, SpatialExtent(cx-500, cx+500, cy-500, cy+500))

fig = Figure(size=(600, 600))
ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true,
           xlabel="x (µm)", ylabel="y (µm)")

# Each call applies the crop independently; no intermediate copies
heatmap!(ax, roi.images["morphology_focus"]; channel=1, colormap=:grays)
poly!(ax,    roi.shapes["cell_boundaries"];  color=:transparent,
                                             strokecolor=:cyan, strokewidth=0.4)
scatter!(ax, roi.points["transcripts"];      color=(:red, 0.25), markersize=1)
fig
```

When a `SpatialElementView{SpatialImage}` is passed to `heatmap!`, the axis
limits are automatically constrained to the view extent — no explicit `limits!`
call is needed.

You can also create element-level views directly:

```julia
pts_view = view(xen.points["transcripts"], extent)
scatter(pts_view; markersize=1)
```

## Multi-channel display

```julia
channel_labels = ["DAPI", "ATP1A1/CD45/E-Cadherin", "18S", "AlphaSMA/Vimentin"]
colormaps      = [:grays, :viridis, :magma, :plasma]

fig = Figure(size=(900, 250))
for (i, (lbl, cmap)) in enumerate(zip(channel_labels, colormaps))
    ax = Axis(fig[1, i]; title=lbl, aspect=DataAspect(),
              yreversed=true, xticksvisible=false, yticksvisible=false)
    heatmap!(ax, img; channel=i, colormap=cmap)
    tightlimits!(ax)
end
fig
```
