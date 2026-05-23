# Quickstart

## Loading data

SpatialOmics uses format tokens to select the reader. Pass the token as the
first argument to `read`:

```julia
using SpatialOmics
import SpatialOmics as SO

# SpatialData OME-Zarr (Xenium, CosMx, Visium, MERFISH — any SpatialData-compatible store)
ds = read(SpatialDataZarr(), "/path/to/experiment.zarr")

# CosMx SMI raw flat-file export
ds = read(CosMx(), "/path/to/cosmx_export/")

# CosMx with morphology images (stitched from per-FOV TIF tiles)
ds = read(CosMx(morphology_dir="/path/to/Morphology2D"), "/path/to/cosmx_export/")
```

## Inspecting a dataset

Use typed accessor functions. They retrieve a named element and verify its type:

```julia
# List all element names
keys(elements(ds))

# Retrieve typed elements
img = images(ds, "morphology_focus")   # SpatialImage
tx  = points(ds, "transcripts")        # SpatialPoints
bnd = shapes(ds, "cell_boundaries")    # SpatialShapes

# Coordinate systems registered in this dataset
coord_systems(ds)

# Gene/feature summary
top_features(tx, 10)          # 10 most frequent genes
count_per_instance(tx)        # transcript count per cell
```

## Defining a region of interest

`SpatialExtent` is a typed bounding box. `view(ds, ext)` returns a
`SpatialDatasetView` — a lazy pointer that applies the spatial filter only when
data are accessed:

```julia
ext = SpatialExtent(4000.0, 5000.0, 1000.0, 2000.0; coord_system="global")
roi = view(ds, ext)

# Accessing elements through the view returns lazy filtered views
tx_roi  = points(roi, "transcripts")     # SpatialElementView{SpatialPoints}
bnd_roi = shapes(roi, "cell_boundaries") # SpatialElementView{SpatialShapes}

# Materialise when you need a concrete copy for computation
tx_mat = collect(tx_roi)   # SpatialPoints with only the filtered rows
```

## Visualization

Load a Makie backend **before** `using SpatialOmics` (or before the first plot):

```julia
using GLMakie      # interactive; or CairoMakie (static), WGLMakie (browser)
```

Standard Makie verbs are extended to accept spatial types directly:

```julia
# Build a composite ROI panel
fig = Figure(size=(600, 600))
ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true,
           xlabel="x (µm)", ylabel="y (µm)")

heatmap!(ax, images(roi, "morphology_focus"); channel=1, colormap=:grays)
poly!(ax,    shapes(roi, "cell_boundaries"); color=:transparent, strokecolor=:cyan, strokewidth=0.4)
scatter!(ax, points(roi, "transcripts");     color=(:red, 0.25), markersize=1)
tightlimits!(ax)
fig
```

See the [Visualization guide](visualization.md) for pyramid-aware images,
multi-channel display, and channel selection.

## Saving

```julia
write!(ds, "/path/to/output.zarr", SpatialDataZarr())
```

`write!` writes to disk and updates the dataset's backing store to the new
location — use this for persistent saves. The resulting Zarr directory is
compatible with Python's SpatialData library.
