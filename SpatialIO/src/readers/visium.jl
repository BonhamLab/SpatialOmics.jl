# SpatialIO/src/readers/visium.jl
# 10x Genomics Visium SpaceRanger output reader.
#
# Expected output directory layout:
#   filtered_feature_bc_matrix.h5
#   spatial/tissue_positions.csv   (or .csv.gz in older SpaceRanger versions)
#   spatial/tissue_hires_image.png
#   spatial/scalefactors_json.json

"""
    VisiumReader

Reader for 10x Genomics Visium SpaceRanger output directories.

Fields
------
- `lazy`          : Bool   — lazy load the HDF5 expression matrix (default true)
- `image_key`     : String — which image resolution to load, "hires" or "lowres"
"""
Base.@kwdef struct VisiumReader <: PlatformReader
    lazy::Bool = true
    image_key::String = "hires"
end

"""
    load_visium(path; lazy=true, image_key="hires") -> SpatialDataset

High-level convenience wrapper around `VisiumReader`.

# Example
```julia
ds = load_visium("/data/spaceranger_output/outs/")
```
"""
function load_visium(path::String; lazy::Bool = true, image_key::String = "hires")
    return load_spatial_data(path, VisiumReader; lazy = lazy, image_key = image_key)
end

function validate_path(::VisiumReader, path::String)
    required = [
        "filtered_feature_bc_matrix.h5",
        joinpath("spatial", "tissue_positions.csv"),
        joinpath("spatial", "scalefactors_json.json"),
    ]
    return all(f -> isfile(joinpath(path, f)), required)
end

function read_data(reader::VisiumReader, path::String)::SpatialDataset
    # TODO Phase 1: implement full Visium reader
    # 1. Read scalefactors_json.json for scale factors
    # 2. Read tissue_positions.csv -> SpatialPoints (spot barcodes + coordinates)
    # 3. Read filtered_feature_bc_matrix.h5 -> SpatialTable (expression matrix)
    # 4. Read tissue_hires_image.png -> SpatialImage
    error("VisiumReader.read_data: not yet implemented (Phase 1 deliverable)")
end
