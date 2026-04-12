# SpatialIO/src/formats/parquet.jl
# Low-level Parquet2.jl helpers for transcript and cell tables.

"""
    read(::TranscriptsParquet, path::String; columns=nothing) -> DataFrame

Read a Xenium/MERFISH transcripts Parquet file. `columns` restricts which
columns are loaded (pass `nothing` to load all).
"""
function read(::TranscriptsParquet, path::String; columns = nothing)::DataFrame
    # TODO Phase 1: implement with Parquet2.jl
    error("read TranscriptsParquet: not yet implemented (Phase 1 deliverable)")
end

"""
    read(::CellsParquet, path::String) -> DataFrame

Read a Xenium cells Parquet file into a DataFrame.
"""
function read(::CellsParquet, path::String)::DataFrame
    # TODO Phase 1: implement with Parquet2.jl
    error("read CellsParquet: not yet implemented (Phase 1 deliverable)")
end
