# SpatialOmics/src/utils/metadata.jl
# Convenience helpers for reading and writing metadata on any SpatialElement
# or SpatialDataset.

get_metadata(x::SpatialDataset) = x.metadata
get_metadata(x::SpatialElement) = x.metadata

function set_metadata!(x::SpatialDataset, key::String, value)
    x.metadata[key] = value
    return x
end
