# SpatialIntegration/src/harmonization/harmonize.jl
# Top-level multi-dataset harmonization.

"""
    harmonize_datasets(datasets::Vector{SpatialDataset};
                       batch_correction::Symbol = :harmony,
                       common_features::Bool = true) -> SpatialDataset

Merge and harmonise a collection of `SpatialDataset` objects into a single
unified dataset, optionally restricting to the intersection of features and
applying batch-effect correction.

Steps performed (planned):
1. Feature standardization — unify gene names, handle synonyms
2. Coordinate system alignment — register datasets to a common frame
3. Batch effect correction — remove technical variation using `batch_correction`
4. Table concatenation — merge expression matrices and observation metadata

Supported `batch_correction` values: `:harmony`, `:scvi`, `:combat`, `:none`.
"""
function harmonize_datasets(
    datasets::Vector{SpatialDataset};
    batch_correction::Symbol = :harmony,
    common_features::Bool = true,
)
    # TODO Phase 3: implement
    error("harmonize_datasets: not yet implemented (Phase 3 deliverable)")
end
