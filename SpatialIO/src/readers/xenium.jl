# SpatialIO/src/readers/xenium.jl
# 10x Genomics Xenium in-situ platform reader.
#
# Expected output directory layout:
#   transcripts.parquet
#   cells.parquet
#   cell_feature_matrix.h5
#   morphology.ome.tiff          (optional)
#   experiment.xenium            (JSON metadata)

"""
    XeniumReader

Reader for 10x Genomics Xenium output directories.

Fields
------
- `lazy`       : Bool — if true, large arrays are opened as DiskArrays (default true)
- `chunk_size` : Int  — target chunk size in rows for Parquet reads (default 100_000)
"""
Base.@kwdef struct XeniumReader <: PlatformReader
    lazy::Bool = true
    chunk_size::Int = 100_000
end

"""
    load_xenium(path; lazy=true) -> SpatialDataset

High-level convenience wrapper around `XeniumReader`.

# Example
```julia
ds = load_xenium("/data/xenium_output/")
```
"""
function load_xenium(path::String; lazy::Bool = true)
    return load_spatial_data(path, XeniumReader; lazy = lazy, validate = true)
end

function validate_path(::XeniumReader, path::String)
    required = ["transcripts.parquet", "cells.parquet", "experiment.xenium"]
    return all(f -> isfile(joinpath(path, f)), required)
end

function read_data(reader::XeniumReader, path::String)::SpatialDataset
    # TODO Phase 1: implement full Xenium reader
    # 1. Parse experiment.xenium (JSON) for metadata
    # 2. Read transcripts.parquet -> SpatialPoints
    #    Xenium columns → canonical: feature_name→gene, cell_id (already),
    #    overlaps_nucleus→cell_compartment ("Nuclear"/"Cytoplasm"),
    #    qv→quality_value; coords from x_location/y_location
    # 3. Read cells.parquet -> SpatialTable (obs)
    # 4. Read cell_feature_matrix.h5 -> SpatialTable (expression matrix)
    # 5. Optionally read morphology.ome.tiff -> SpatialImage
    error("XeniumReader.read_data: not yet implemented (Phase 1 deliverable)")
end
