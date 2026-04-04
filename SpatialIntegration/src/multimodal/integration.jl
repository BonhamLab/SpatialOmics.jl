# SpatialIntegration/src/multimodal/integration.jl
# Multi-modal integration (e.g. spatial + single-cell, spatial + proteomics).

"""
    integrate_modalities(datasets::Vector{SpatialDataset};
                         method::Symbol = :mofa,
                         n_factors::Int = 10) -> SpatialDataset

Integrate multiple omics modalities present across the input `datasets` into
a joint low-dimensional embedding stored in the returned dataset.

Supported `method` values (planned):
- `:mofa`   — MOFA+ multi-omics factor analysis (via PythonCall)
- `:mcia`   — Multiple Co-Inertia Analysis
- `:concat` — simple feature concatenation (baseline)
"""
function integrate_modalities(
    datasets::Vector{SpatialDataset};
    method::Symbol = :mofa,
    n_factors::Int = 10,
)
    # TODO Phase 4: implement
    error("integrate_modalities: not yet implemented (Phase 4 deliverable)")
end

"""
    joint_embedding(datasets::Vector{SpatialDataset};
                    method::Symbol = :umap,
                    n_components::Int = 2) -> Matrix{Float32}

Compute a joint low-dimensional embedding across all datasets and modalities.
Returns an N×`n_components` matrix where N is the total number of observations.
"""
function joint_embedding(
    datasets::Vector{SpatialDataset};
    method::Symbol = :umap,
    n_components::Int = 2,
)
    # TODO Phase 4: implement
    error("joint_embedding: not yet implemented (Phase 4 deliverable)")
end
