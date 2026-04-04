# SpatialIntegration/src/reference/reference_mapping.jl
# Map a query SpatialDataset onto an annotated reference atlas.

"""
    integrate_reference(query::SpatialDataset,
                        reference::SpatialDataset;
                        method::Symbol = :tangram,
                        transfer_obs_keys::Vector{String} = ["cell_type"]) -> SpatialDataset

Map cells/spots in `query` to the most similar cells in `reference` using
the selected `method`, and optionally transfer observation-level annotations
listed in `transfer_obs_keys`.

Returns a new `SpatialDataset` equivalent to `query` augmented with the
transferred annotations.

Supported `method` values (planned):
- `:tangram`   — Tangram spatial mapping (via PythonCall)
- `:seurat_v3` — Seurat v3 anchor-based integration (via RCall)
- `:scvi`      — scVI reference mapping
"""
function integrate_reference(
    query::SpatialDataset,
    reference::SpatialDataset;
    method::Symbol = :tangram,
    transfer_obs_keys::Vector{String} = ["cell_type"],
)
    # TODO Phase 3: implement
    error("integrate_reference: not yet implemented (Phase 3 deliverable)")
end
