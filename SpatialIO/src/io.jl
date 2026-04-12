# SpatialIO/src/io.jl
#
# Unified read/write dispatch API.
#
# All format-level I/O goes through two unexported functions:
#
#   SpatialIO.read(format, path; kwargs...)
#   SpatialIO.write(data, path, format; kwargs...)
#
# where `format` is an instance of one of the StorageFormat singleton types
# defined below.  Import the module qualified to use them:
#
#   import SpatialIO as SIO
#   ds = SIO.read(SIO.Zarr(), "/data/experiment.zarr")
#   SIO.write(ds, "/out/experiment.h5", SIO.NativeH5())

# ─── Format type hierarchy ────────────────────────────────────────────────────

"""
Abstract supertype for all on-disk format tokens.
Subtypes are singleton structs used as the first (read) or third (write)
argument to `SpatialIO.read` / `SpatialIO.write` to select the format.
"""
abstract type StorageFormat end

"""
    Zarr()

SpatialData-compatible OME-ZARR store (Zarr v3).
Backed by the native `from_spatialdata` / `to_spatialdata` implementation.
"""
struct Zarr <: StorageFormat end

"""
    NativeH5()

Native SpatialIO HDF5 serialisation (round-trips all element types).
Not the same as the 10x Genomics cell-feature-barcode HDF5 format.
"""
struct NativeH5 <: StorageFormat end

"""
    AnnData()

AnnData-compatible `.h5ad` export (write-only).
"""
struct AnnData <: StorageFormat end

"""
    GeoJSON()

GeoJSON FeatureCollection export for `SpatialShapes` (write-only).
"""
struct GeoJSON <: StorageFormat end

# ── Vendor / platform file formats ───────────────────────────────────────────

"""
    TranscriptsParquet()

Xenium or MERFISH per-transcript Parquet file.
`read` returns a `DataFrame`.
"""
struct TranscriptsParquet <: StorageFormat end

"""
    CellsParquet()

Xenium cells-summary Parquet file.
`read` returns a `DataFrame`.
"""
struct CellsParquet <: StorageFormat end

"""
    CellFeatureMatrixH5()

10x Genomics `filtered_feature_bc_matrix.h5` / `cell_feature_matrix.h5`.
`read` returns `(matrix::AbstractMatrix, barcodes::Vector{String}, features::DataFrame)`.
"""
struct CellFeatureMatrixH5 <: StorageFormat end

"""
    TissuePositionsCSV()

Visium `tissue_positions.csv` (or `.csv.gz`) file.
`read` returns a `DataFrame`.
"""
struct TissuePositionsCSV <: StorageFormat end

# ─── read ─────────────────────────────────────────────────────────────────────

"""
    read(path::String; kwargs...) -> SpatialDataset

Read a spatial dataset from `path`, guessing the format from the file extension
or directory structure:
- `.zarr` / directory containing `zarr.json` → `Zarr()`
- `.h5` / `.hdf5`                            → `NativeH5()`

Keyword arguments are forwarded to the selected format's `read` method.
"""
function read(path::String; kwargs...)
    if endswith(path, ".zarr") || isdir(path) && isfile(joinpath(path, "zarr.json"))
        return read(Zarr(), path; kwargs...)
    elseif endswith(path, ".h5") || endswith(path, ".hdf5")
        return read(NativeH5(), path; kwargs...)
    else
        error("read: cannot guess format for \"$path\". Pass an explicit format token, e.g. read(Zarr(), path).")
    end
end

"""
    read(::Zarr, path::String; use_python::Bool=false) -> SpatialDataset

Read a SpatialData OME-ZARR store.

# Example
```julia
import SpatialIO as SIO
ds = SIO.read(SIO.Zarr(), "/data/experiment.zarr")
```
"""
function read(::Zarr, path::String; kwargs...)::SpatialDataset
    from_spatialdata(path; kwargs...)
end

"""
    read(::NativeH5, path::String) -> SpatialDataset

Read a `SpatialDataset` from a native SpatialIO HDF5 file.

# Example
```julia
import SpatialIO as SIO
ds = SIO.read(SIO.NativeH5(), "my_dataset.h5")
```
"""
function read(::NativeH5, path::String)::SpatialDataset
    read_hdf5(path)
end

# Format-specific read methods are defined in their respective format files:
#   formats/parquet.jl  — TranscriptsParquet, CellsParquet
#   formats/hdf5.jl     — CellFeatureMatrixH5
#   formats/csv.jl      — TissuePositionsCSV

# ─── write ────────────────────────────────────────────────────────────────────

"""
    write(ds::SpatialDataset, path::String, ::Zarr; use_python::Bool=false)

Write `ds` to a SpatialData OME-ZARR store.

# Example
```julia
import SpatialIO as SIO
SIO.write(ds, "/out/experiment.zarr", SIO.Zarr())
```
"""
function write(ds::SpatialDataset, path::String, ::Zarr; kwargs...)
    to_spatialdata(ds, path; kwargs...)
end

"""
    write(ds::SpatialDataset, path::String, ::NativeH5)

Serialise `ds` to a native SpatialIO HDF5 file.

# Example
```julia
import SpatialIO as SIO
SIO.write(ds, "my_dataset.h5", SIO.NativeH5())
```
"""
function write(ds::SpatialDataset, path::String, ::NativeH5)
    write_hdf5(ds, path)
end

# Format-specific write methods are defined in their respective format files:
#   formats/anndata.jl  — AnnData
#   formats/geojson.jl  — GeoJSON
