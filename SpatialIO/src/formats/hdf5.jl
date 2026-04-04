# SpatialIO/src/formats/hdf5.jl
# Helpers for reading 10x HDF5 feature-barcode matrices.

"""
    read_cell_feature_matrix_h5(path::String; lazy::Bool=true)
        -> (matrix::AbstractMatrix, barcodes::Vector{String}, features::DataFrame)

Read a `filtered_feature_bc_matrix.h5` or `cell_feature_matrix.h5` file.
Returns the sparse count matrix, barcode strings, and a feature DataFrame.
When `lazy=true`, the matrix is memory-mapped rather than fully loaded.
"""
function read_cell_feature_matrix_h5(path::String; lazy::Bool = true)
    # TODO Phase 1: implement with HDF5.jl
    # Standard 10x HDF5 layout:
    #   /matrix/barcodes        (string dataset)
    #   /matrix/features/id     (string dataset)
    #   /matrix/features/name   (string dataset)
    #   /matrix/data            (int32 — CSC sparse values)
    #   /matrix/indices         (int32 — CSC row indices)
    #   /matrix/indptr          (int32 — CSC column pointers)
    #   /matrix/shape           (int32[2])
    error("read_cell_feature_matrix_h5: not yet implemented (Phase 1 deliverable)")
end
