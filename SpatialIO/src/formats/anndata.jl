# SpatialIO/src/formats/anndata.jl
# Export a SpatialDataset (or a single SpatialTable) to AnnData .h5ad format.

"""
    write(ds::SpatialDataset, path::String, ::AnnData; table_key::String="table")

Export the named SpatialTable from `ds` as an AnnData-compatible HDF5 file
at `path`. Spatial coordinates are stored in `obsm["spatial"]`.
"""
function write(ds::SpatialDataset, path::String, ::AnnData; table_key::String = "table")
    # TODO Phase 1: implement
    # AnnData HDF5 layout:
    #   /X           — count matrix (dense or sparse)
    #   /obs         — per-observation annotations
    #   /var         — per-variable (gene) annotations
    #   /obsm/spatial — N×2 spatial coordinate matrix
    error("write AnnData: not yet implemented (Phase 1 deliverable)")
end
