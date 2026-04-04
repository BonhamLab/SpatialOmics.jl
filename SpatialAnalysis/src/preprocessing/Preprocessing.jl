"""
    Preprocessing

Count normalisation, batch effect correction, and missing-value imputation
for spatial transcriptomics datasets. Functions are designed to be pure and
composable, operating on `SpatialDataset` objects.
"""
module Preprocessing

using DataFrames
using Statistics
using SpatialOmicsBase

export normalize_counts, correct_batch_effects, impute_missing

"""
    normalize_counts(ds::SpatialDataset;
                     method::Symbol = :scran,
                     target_sum::Real = 1e4) -> SpatialDataset

Normalise raw count data in `ds`.

Supported `method` values (planned):
- `:scran`    — pooling-based size-factor normalisation (via PythonCall bridge)
- `:total`    — divide each cell by its total count, multiply by `target_sum`
- `:log1p`    — log1p-transform after total normalisation
"""
function normalize_counts(
    ds::SpatialDataset;
    method::Symbol = :scran,
    target_sum::Real = 1e4,
)
    # TODO Phase 2: implement
    error("normalize_counts: not yet implemented (Phase 2 deliverable)")
end

"""
    correct_batch_effects(ds::SpatialDataset;
                          batch_key::String = "sample",
                          method::Symbol = :harmony) -> SpatialDataset

Remove technical batch variation across slides or samples.

Supported `method` values (planned):
- `:harmony`  — Harmony integration (via PythonCall)
- `:combat`   — ComBat parametric batch correction
- `:scvi`     — scVI deep generative model (via PythonCall)
"""
function correct_batch_effects(
    ds::SpatialDataset;
    batch_key::String = "sample",
    method::Symbol = :harmony,
)
    # TODO Phase 2: implement
    error("correct_batch_effects: not yet implemented (Phase 2 deliverable)")
end

"""
    impute_missing(ds::SpatialDataset;
                   method::Symbol = :tangram,
                   reference::Union{Nothing, SpatialDataset} = nothing)
        -> SpatialDataset

Impute missing gene expression values using spatial context.

Supported `method` values (planned):
- `:tangram`  — Tangram transcript-level mapping (requires reference)
- `:gimvi`    — GIMVI spatial imputation
- `:knn`      — simple k-NN smoothing over spatial neighbours
"""
function impute_missing(
    ds::SpatialDataset;
    method::Symbol = :tangram,
    reference::Union{Nothing, SpatialDataset} = nothing,
)
    # TODO Phase 2: implement
    error("impute_missing: not yet implemented (Phase 2 deliverable)")
end

end # module Preprocessing
