# SpatialIntegration/src/reference/deconvolution.jl
# Deconvolute bulk-like spots (Visium) into cell-type proportions.

"""
    deconvolute_spots(query::SpatialDataset,
                      reference::SpatialDataset;
                      method::Symbol = :rctd,
                      cell_type_key::String = "cell_type") -> DataFrame

Estimate the fractional composition of cell types within each spatial
observation (spot) of `query` using single-cell `reference` data.

Returns a DataFrame with one row per spot and one column per cell type,
containing the estimated cell-type proportions (values sum to ~1 per row).

Supported `method` values (planned): `:rctd`, `:card`, `:spotlight`.
"""
function deconvolute_spots(
    query::SpatialDataset,
    reference::SpatialDataset;
    method::Symbol = :rctd,
    cell_type_key::String = "cell_type",
)
    # TODO Phase 3: implement
    error("deconvolute_spots: not yet implemented (Phase 3 deliverable)")
end
