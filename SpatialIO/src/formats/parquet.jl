# SpatialIO/src/formats/parquet.jl
# Low-level Parquet2.jl helpers for transcript and cell tables.

"""
    read_transcripts_parquet(path::String; columns=nothing) -> DataFrame

Read a Xenium/MERFISH transcripts Parquet file. `columns` restricts which
columns are loaded (pass `nothing` to load all).
"""
function read_transcripts_parquet(path::String; columns = nothing)
    # TODO Phase 1: implement with Parquet2.jl
    error("read_transcripts_parquet: not yet implemented (Phase 1 deliverable)")
end

"""
    read_cells_parquet(path::String) -> DataFrame

Read a Xenium cells Parquet file into a DataFrame.
"""
function read_cells_parquet(path::String)
    # TODO Phase 1: implement with Parquet2.jl
    error("read_cells_parquet: not yet implemented (Phase 1 deliverable)")
end
