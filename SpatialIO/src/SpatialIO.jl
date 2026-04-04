"""
    SpatialIO.jl

Platform I/O layer for the STX_DEV spatial transcriptomics framework.
Provides pluggable `PlatformReader` subtypes for Xenium, Visium, CosMx, and
MERFISH vendor formats, plus round-trip interoperability with the Python
SpatialData ecosystem via native Zarr/HDF5/Arrow/Parquet
backends. All readers return `SpatialDataset` objects from SpatialOmicsBase.jl.

## Format-level I/O

Use the unexported `read` / `write` functions with format type tokens:

```julia
import SpatialIO as SIO

ds = SIO.read(SIO.Zarr(),     "/data/experiment.zarr")
ds = SIO.read(SIO.NativeH5(), "experiment.h5")

SIO.write(ds, "/out/experiment.zarr", SIO.Zarr())
SIO.write(ds, "experiment.h5",        SIO.NativeH5())
SIO.write(ds, "table.h5ad",           SIO.AnnData())
```
"""
module SpatialIO

using HDF5
using Arrow
using CSV
using Parquet2
using JSON3
using Printf: @sprintf
using CodecZstd
using DiskArrays
using SparseArrays
using DataFrames
using TiffImages
using SpatialOmicsBase

# ---------------------------------------------------------------------------
# Nothing is exported from SpatialIO.
#
# All public names are surfaced through the SpatialOmics umbrella:
#   import SpatialOmics as IO
#   IO.load_cosmx(path)
#   IO.read(IO.Zarr(), path)
#   IO.write(ds, path, IO.NativeH5())
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Includes
# ---------------------------------------------------------------------------

include("readers/platform_reader.jl")
include("readers/xenium.jl")
include("readers/visium.jl")
include("readers/cosmx.jl")
include("readers/merfish.jl")
include("spatialdata/zarr3.jl")
include("spatialdata/from_spatialdata.jl")
include("spatialdata/to_spatialdata.jl")
include("formats/parquet.jl")
include("formats/hdf5.jl")
include("formats/csv.jl")
include("formats/anndata.jl")
include("formats/geojson.jl")
include("hdf5_backend.jl")
include("zarr_backend.jl")
include("io.jl")

end # module SpatialIO
