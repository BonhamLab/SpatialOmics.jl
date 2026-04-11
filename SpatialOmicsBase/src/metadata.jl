# SpatialOmics/src/utils/metadata.jl
# Convenience helpers for reading and writing metadata on any SpatialElement
# or SpatialDataset.

# TODO:: This seems like dead code - pay attention to memories,
# get_* and set_* method definitions are anti-patterns.
# I think there other accessors defined elsewhere now.
get_metadata(x::SpatialDataset) = x.metadata
get_metadata(x::SpatialElement) = x.metadata

function set_metadata!(x::SpatialDataset, key::String, value)
    x.metadata[key] = value
    return x
end
