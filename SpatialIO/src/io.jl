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

"""
    read(::TranscriptsParquet, path::String; columns=nothing) -> DataFrame

Read a Xenium / MERFISH per-transcript Parquet file.
Pass `columns` to restrict which columns are loaded.
"""
function read(::TranscriptsParquet, path::String; columns = nothing)::DataFrame
    read_transcripts_parquet(path; columns)
end

"""
    read(::CellsParquet, path::String) -> DataFrame

Read a Xenium cells-summary Parquet file into a DataFrame.
"""
function read(::CellsParquet, path::String)::DataFrame
    read_cells_parquet(path)
end

"""
    read(::CellFeatureMatrixH5, path::String; lazy::Bool=true)
        -> (matrix::AbstractMatrix, barcodes::Vector{String}, features::DataFrame)

Read a 10x Genomics cell-feature-barcode HDF5 file.
"""
function read(::CellFeatureMatrixH5, path::String; lazy::Bool = true)
    read_cell_feature_matrix_h5(path; lazy)
end

"""
    read(::TissuePositionsCSV, path::String) -> DataFrame

Parse a Visium `tissue_positions.csv` (or `.csv.gz`) file.
"""
function read(::TissuePositionsCSV, path::String)::DataFrame
    read_tissue_positions_csv(path)
end

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

"""
    write(ds::SpatialDataset, path::String, ::AnnData; table_key::String="table")

Export the named `SpatialTable` from `ds` as an AnnData `.h5ad` file.
"""
function write(ds::SpatialDataset, path::String, ::AnnData; table_key::String = "table")
    write_anndata(ds, path; table_key)
end

"""
    write(shapes::SpatialShapes, path::String, ::GeoJSON)

Write `shapes` to a GeoJSON FeatureCollection file.
"""
function write(shapes::SpatialShapes, path::String, ::GeoJSON)
    write_geojson(shapes, path)
end
