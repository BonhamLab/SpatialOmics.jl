# I/O

SpatialOmics uses format tokens to select import and snapshot formats. Native
dataset changes are persisted with [`save!`](@ref).

## Format tokens

```@docs
SpatialDataZarr
CosMx
native_store_version
```

## Reading

`Base.read` is extended for spatial format tokens:

```julia
# Native SpatialOmics Zarr, or a supported Python-written SpatialData store
ds = read(SpatialDataZarr(), "/path/to/experiment.zarr")

# CosMx SMI raw flat-file export
ds = read(CosMx(), "/path/to/cosmx_export/")

# CosMx with tissue images
ds = read(CosMx(morphology_dir="/path/to/Morphology2D"), "/path/to/export/")
```

Native stores carry an explicit schema version. Stores written before versioning,
or with an unsupported version, are rejected with a rebuild instruction; reads never
silently migrate or reinterpret them. Use `native_store_version(path)` to inspect a
store before opening it. A future migration API can upgrade a store only when called
explicitly.

## Writing

```@docs
write!
```

`write!` is retained as a compatibility spelling for a complete native save
and rebind. `Base.write` writes a snapshot without updating the active backing
location or clearing its dirty state.
