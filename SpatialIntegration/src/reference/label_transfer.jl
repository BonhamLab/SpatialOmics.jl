# SpatialIntegration/src/reference/label_transfer.jl
# Transfer discrete cell-type labels from a reference to a query dataset.

"""
    transfer_labels(query::SpatialDataset,
                    reference::SpatialDataset,
                    label_key::String = "cell_type";
                    method::Symbol = :knn,
                    k::Int = 15) -> Vector{String}

Predict cell-type labels for each observation in `query` based on similarity
to labelled observations in `reference`. Returns a vector of predicted label
strings with one entry per observation in `query`.

Supported `method` values (planned): `:knn`, `:svm`, `:tangram`.
"""
function transfer_labels(
    query::SpatialDataset,
    reference::SpatialDataset,
    label_key::String = "cell_type";
    method::Symbol = :knn,
    k::Int = 15,
)
    # TODO Phase 3: implement
    error("transfer_labels: not yet implemented (Phase 3 deliverable)")
end
