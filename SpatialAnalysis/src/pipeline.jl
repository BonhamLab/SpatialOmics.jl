# SpatialAnalysis/src/pipeline.jl
# High-level pipeline helpers that chain QC and preprocessing steps.

"""
    default_qc_pipeline(; min_transcripts=10, min_genes=5) -> Vector{Function}

Return a default ordered list of QC steps suitable for most Xenium/Visium
datasets. Steps are plain functions that each accept and return a
`SpatialDataset`.
"""
function default_qc_pipeline(; min_transcripts::Int = 10, min_genes::Int = 5)
    return [
        ds -> filter_spots(ds; min_transcripts = min_transcripts, min_genes = min_genes),
        ds -> normalize_counts(ds; method = :total),
    ]
end

"""
    run_qc_pipeline(ds::SpatialDataset,
                    steps::Vector{Function} = default_qc_pipeline()) -> SpatialDataset

Thread `ds` through each step function in `steps` sequentially, returning
the final processed `SpatialDataset`.
"""
function run_qc_pipeline(
    ds::SpatialDataset,
    steps::Vector{Function} = default_qc_pipeline(),
)
    result = ds
    for step in steps
        result = step(result)
    end
    return result
end
