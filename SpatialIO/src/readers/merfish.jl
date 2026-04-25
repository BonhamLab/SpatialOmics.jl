# SpatialIO/src/readers/merfish.jl
# Vizgen MERSCOPE / MERFISH output reader.
#
# Key files: detected_transcripts.csv, cell_metadata.csv,
#            cell_by_gene.csv, images/mosaic_*.tif, cell_boundaries.parquet

"""
    MerfishReader

Reader for Vizgen MERSCOPE (MERFISH) output directories.

Fields
------
- `lazy`      : Bool   — lazy image loading (default true)
- `z_slice`   : Int    — which z-slice to load for 3-D datasets (default 2,
                         typically the mid-z plane)
"""
Base.@kwdef struct MerfishReader <: PlatformReader
    lazy::Bool = true
    z_slice::Int = 2
end

"""
    load_merfish(path; lazy=true, z_slice=2) -> SpatialDataset

High-level convenience wrapper around `MerfishReader`.
"""
function load_merfish(path::String; lazy::Bool = true, z_slice::Int = 2)
    return load_spatial_data(path, MerfishReader; lazy = lazy, z_slice = z_slice)
end

function read_data(reader::MerfishReader, path::String)::SpatialDataset
    # TODO Phase 1: implement MERFISH reader
    # 1. Read detected_transcripts.csv -> SpatialPoints
    #    MERFISH columns → canonical: gene (already), cell_id (already),
    #    global_x/global_y → coords; barcode_id and other platform cols kept as-is
    # 2. Read cell_metadata.csv + cell_by_gene.csv -> SpatialTable
    # 3. Read mosaic TIFF images -> SpatialImage (selected z-slice)
    # 4. Read cell_boundaries.parquet -> SpatialShapes
    error("MerfishReader.read_data: not yet implemented (Phase 1 deliverable)")
end
