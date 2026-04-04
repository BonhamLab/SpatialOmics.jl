# SpatialIntegration/src/harmonization/feature_standardization.jl
# Standardize gene/feature names across datasets from different platforms.

"""
    standardize_features(datasets::Vector{SpatialDataset};
                         id_type::Symbol = :symbol,
                         species::String = "human") -> Vector{SpatialDataset}

Normalize feature (gene) identifiers across datasets to a common namespace.
Handles Ensembl ID ↔ HGNC symbol conversion and synonym resolution.

`id_type` options: `:symbol` (HGNC gene symbol), `:ensembl` (Ensembl gene ID).
"""
function standardize_features(
    datasets::Vector{SpatialDataset};
    id_type::Symbol = :symbol,
    species::String = "human",
)
    # TODO Phase 3: implement (may use BioMart via PythonCall or a lookup table)
    error("standardize_features: not yet implemented (Phase 3 deliverable)")
end
