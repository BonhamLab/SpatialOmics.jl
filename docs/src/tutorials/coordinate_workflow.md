# Place an FOV in slide coordinates

```@setup coordinate-workflow
using CairoMakie
using Markdown
CairoMakie.activate!(type="svg")
set_theme!(Theme(
    fontsize=15,
    Figure=(; backgroundcolor=:white),
    Axis=(; xgridvisible=false, ygridvisible=false),
))
```

Imaging assays often report transcript positions in FOV-local pixels while
slide overlays use a physical global coordinate system. SpatialOmics keeps
those spaces named and resolves the transformation path explicitly.

## Define the coordinate spaces

This example has 0.5 µm pixels. The FOV origin lies at `(1000, 250)` µm in the
slide coordinate system.

```@example coordinate-workflow
using SpatialOmics
import SpatialOmics as SO

dataset = SpatialDataset()
push!(dataset, CoordinateSystem("fov_1_px"; units=("px", "px")))
push!(dataset, CoordinateSystem("fov_1_um"; units=("µm", "µm")))
push!(dataset, CoordinateSystem("slide_um"; units=("µm", "µm")))

pixel_size = SO.scaling(0.5, 0.5, "fov_1_px", "fov_1_um")
placement = SO.translation(1000.0, 250.0, "fov_1_um", "slide_um")
push!(dataset, pixel_size)
push!(dataset, placement)

coord_systems(dataset)
```

## Resolve and apply the path

`transform` resolves a path through the named graph. `apply` returns a new
element and leaves the local data unchanged.

```@example coordinate-workflow
local_transcripts = SpatialPoints(
    [Point2f(0, 0), Point2f(20, 40), Point2f(60, 20),
     Point2f(80, 70), Point2f(100, 50)];
    coord_system="fov_1_px",
)

to_slide = transform(dataset, "fov_1_px", "slide_um")
slide_transcripts = apply(to_slide, local_transcripts)

(
    local_coordinates=coords(local_transcripts),
    slide=coords(slide_transcripts),
    destination=coord_system(slide_transcripts),
)
```

The same observations retain their identity as the FOV is scaled from pixels
to micrometres and translated into its slide position.

```@eval coordinate-workflow
point_colors = Makie.wong_colors()[1:length(local_transcripts)]
figure = Figure(size=(820, 400))
local_axis = Axis(figure[1, 1]; title="FOV-local pixels",
                  xlabel="x (px)", ylabel="y (px)")
slide_axis = Axis(figure[1, 2]; title="Placed on slide",
                  xlabel="x (µm)", ylabel="y (µm)")
scatter!(local_axis, local_transcripts; color=point_colors, markersize=14)
scatter!(slide_axis, slide_transcripts; color=point_colors, markersize=14)
poly!(local_axis, [Rect2f(0, 0, 110, 80)]; color=(:steelblue, 0.08),
      strokecolor=:steelblue, strokewidth=2)
poly!(slide_axis, [Rect2f(1000, 250, 55, 40)]; color=(:steelblue, 0.08),
      strokecolor=:steelblue, strokewidth=2)
Label(figure[2, 1:2], "scale × 0.5, then translate + (1000, 250)",
      tellheight=true, justification=:center)
xlims!(local_axis, -8, 118); ylims!(local_axis, -8, 88)
xlims!(slide_axis, 996, 1059); ylims!(slide_axis, 246, 294)
save("coordinate-workflow.svg", figure)
Markdown.parse("![FOV coordinates before and after transformation](coordinate-workflow.svg)")
```

Affine edges can also be resolved in reverse. Unsupported paths fail instead
of silently combining incompatible coordinates.

```@example coordinate-workflow
to_pixels = transform(dataset, "slide_um", "fov_1_px")
apply(to_pixels, coords(slide_transcripts))
```

For attached elements, prefer the pure `apply` form unless the dataset itself
should change. Package-mediated mutation with `apply!` marks an attached
element dirty; other direct field mutation must be wrapped in `edit!` or
followed by `touch!`.

```@example coordinate-workflow
close(dataset; discard=true) # hide
nothing # hide
```
