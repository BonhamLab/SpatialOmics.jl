# I/O

SpatialOmics uses format tokens to select the read/write backend. Pass the
token as the first argument to `read` or `write!`.

## Format tokens

```@docs
SpatialDataZarr
CosMx
```

## Reading

`Base.read` is extended for spatial format tokens:

```julia
# SpatialData OME-Zarr (auto-detects Julia vs Python-written stores)
ds = read(SpatialDataZarr(), "/path/to/experiment.zarr")

# CosMx SMI raw flat-file export
ds = read(CosMx(), "/path/to/cosmx_export/")

# CosMx with tissue images
ds = read(CosMx(morphology_dir="/path/to/Morphology2D"), "/path/to/export/")
```

## Writing

```@docs
write!
```

`Base.write` (without `!`) is also defined and writes to disk without updating
the dataset's backing store location. Prefer `write!` for persistent saves.
