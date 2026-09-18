# CosMx workflow

CosMx SMI exports a flat-file directory with per-FOV transcripts, cell
segmentation polygons, and optional tissue images (Morphology2D TIF tiles).

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
bnd = shapes(ds, "cell_boundaries")

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
coord_systems(ds)   # ["fov_1", "fov_2", ..., "global"]

# Resolve a transform from a FOV to global space
t = transform(ds, "fov_1", "global")

# Apply to transform an element between spaces
tx_global = apply(t, points(ds, "transcripts_fov_1"))
```

## Spatial filtering

Use `SpatialExtent` or `SpatialROI` to define a region of interest. Views
are lazy — no data is copied:

```julia
ext = SpatialExtent(5000.0, 7000.0, 3000.0, 5000.0; coord_system="global")
roi = view(ds, ext)

# Filter transcripts and shapes to the ROI
tx_roi  = points(roi, "transcripts")
bnd_roi = shapes(roi, "cell_boundaries")

# Subsampled scatter for quick overview
scatter!(ax, subsample(collect(tx_roi), 50_000); markersize=1)
```

## Visualisation

```julia
using CairoMakie

ext = SpatialExtent(5000.0, 6000.0, 3000.0, 4000.0; coord_system="global")
roi = view(ds, ext)

fig = Figure(size=(600, 600))
ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true)

# Tissue image — rescaled for display
image!(ax, scaleminmax(channel(images(roi, "morphology"), 1)))
# Cell boundaries
poly!(ax, shapes(roi, "cell_boundaries"); color=:transparent, strokecolor=:cyan, strokewidth=0.3)
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
