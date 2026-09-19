# Place an FOV in slide coordinates

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
    [Point2f(0, 0), Point2f(20, 40)];
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
