# CosMx workflow

CosMx SMI exports a flat-file directory with per-FOV transcripts, cell
segmentation polygons, and optional tissue images (Morphology2D TIF tiles).

For a reproducible public input, Bruker publishes a [CosMx Human Lymph Node
FFPE dataset](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-human-lymph-node-ffpe-dataset/)
with transcript coordinates, cell metadata, FOV positions, polygons, and
images. Full public releases are too large for routine documentation builds;
rendered examples should be generated from a documented subset and committed
in the same way as the Xenium and Visium tutorial figures.

## Loading

```julia
using SpatialOmics

# Basic load — transcripts and cell boundaries only
ds = read(CosMx(), "/path/to/cosmx_export/")

# Include tissue images (stitched from Morphology2D TIF tiles)
ds = read(CosMx(morphology_dir="/path/to/Morphology2D"), "/path/to/cosmx_export/")

# Cache to disk for faster subsequent loads
save!(ds; path="/path/to/cache.zarr")
ds2 = read(SpatialDataZarr(), "/path/to/cache.zarr")
```

## Exploring elements

```julia
keys(elements(ds))   # list all loaded elements

tx  = points(ds, "transcripts")
bnd = shapes(ds, "cells")

# Top expressed genes
top_features(tx, 20)

# Transcripts per cell (excludes unassigned; instance_id == 0)
counts = count_per_instance(tx)
```

## Coordinate systems

Each FOV is registered as a separate `CoordinateSystem`. All FOVs are
registered in a shared global coordinate system, with per-FOV affine
transforms in the dataset's transform graph.

```julia
# List all registered coordinate systems
coord_systems(ds)   # ["global_px", "fov_1_px", "fov_2_px", ...]

# Resolve a transform from a FOV to global space
t = transform(ds, "fov_1_px", "global_px")

# Apply to transform an element between spaces
local_points = SpatialPoints([Point2f(10, 20)]; coord_system="fov_1_px")
tx_global = apply(t, local_points)
```

## Spatial filtering

Use `SpatialExtent` or `SpatialROI` to define a region of interest. Views
are lazy — no data is copied:

```julia
ext = SpatialExtent(5000.0, 7000.0, 3000.0, 5000.0; coord_system="global_px")
roi = view(ds, ext)

# Filter transcripts and shapes to the ROI
tx_roi  = points(roi, "transcripts")
bnd_roi = shapes(roi, "cells")

# Subsampled scatter for quick overview
scatter!(ax, subsample(collect(tx_roi), 50_000); markersize=1)
```

Each CosMx FOV is also registered as an acquisition source. Source views use
the FOV recorded by the instrument rather than footprint geometry:

```julia
sources(ds)                    # ["fov_1_px", "fov_2_px", ...]
fov2 = view(ds, "fov_2_px")
tx_fov2 = points(fov2, "transcripts")
```

If two FOV footprints overlap, `tx_fov2` contains only transcripts acquired in
FOV 2. A user-drawn `SpatialROI` over the same overlap contains transcripts
from both FOVs. The `z` and `CellComp` transcript annotations are available as
`features(tx, :z)` and `features(tx, :CellComp)`.

## Visualisation

```julia
using CairoMakie

ext = SpatialExtent(5000.0, 6000.0, 3000.0, 4000.0; coord_system="global_px")
roi = view(ds, ext)

fig = Figure(size=(600, 600))
ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true)

# Tissue image — rescaled for display
image!(ax, scaleminmax(channel(images(roi, "morphology"), 1)))
# Cell boundaries
poly!(ax, shapes(roi, "cells"); color=:transparent, strokecolor=:cyan, strokewidth=0.3)
# Top gene transcripts
for gene in top_features(points(roi, "transcripts"), 3)
    scatter!(ax, coords(points(roi, "transcripts"), gene); label=gene, markersize=2)
end
axislegend(ax; position=:lt)
tightlimits!(ax)
fig
```

## Native persistence

Save the assembled dataset in the native SpatialOmics Zarr layout for later
Julia workflows. Python handoff requires an explicit SpatialData export path.

```julia
save!(ds; path="/path/to/output.zarr")
```
