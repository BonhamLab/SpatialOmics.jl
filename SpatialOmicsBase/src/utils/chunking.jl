# SpatialOmics/src/utils/chunking.jl
# Helpers for automatic chunk-size selection when wrapping arrays in DiskArrays.

"""
    default_chunk_size(sz::Tuple{Vararg{Int}}, element_bytes::Int=4) -> Tuple

Return a chunk size tuple that targets ~64 MB per chunk given the full array
size `sz` and per-element byte count `element_bytes`.
This is a placeholder; the heuristic will be refined in Phase 1.
"""
function default_chunk_size(sz::Tuple{Vararg{Int}}, element_bytes::Int = 4)
    # TODO Phase 1: implement adaptive chunking based on available memory
    target_bytes = 64 * 1024 * 1024  # 64 MB
    ndims = length(sz)
    chunk_dim = round(Int, (target_bytes / element_bytes)^(1 / ndims))
    return ntuple(i -> min(sz[i], chunk_dim), ndims)
end
