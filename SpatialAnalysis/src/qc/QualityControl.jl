"""
    QualityControl

Quality control functions for spatial transcriptomics datasets.
Provides spot/cell filtering, artifact detection, and spatially-aware
QC metrics that account for tissue architecture.

All functions accept a `SpatialDataset` and return either a filtered
`SpatialDataset` or a QC result DataFrame, making them composable in
analysis pipelines.
"""
module QualityControl

using DataFrames
using Statistics
using SpatialOmicsBase

export filter_spots, detect_artifacts, spatially_aware_qc

"""
    filter_spots(ds::SpatialDataset;
                 min_transcripts::Int = 10,
                 max_transcripts::Int = typemax(Int),
                 min_genes::Int = 5,
                 in_tissue::Bool = true) -> SpatialDataset

Remove low-quality spots/cells from `ds` based on per-spot transcript count
and gene detection thresholds. For Visium data, `in_tissue=true` (default)
drops spots not overlapping tissue.

Returns a new `SpatialDataset` with the filtered observations; the original
is unchanged.
"""
function filter_spots(
    ds::SpatialDataset;
    min_transcripts::Int = 10,
    max_transcripts::Int = typemax(Int),
    min_genes::Int = 5,
    in_tissue::Bool = true,
)
    # TODO Phase 2: implement
    error("filter_spots: not yet implemented (Phase 2 deliverable)")
end

"""
    detect_artifacts(ds::SpatialDataset;
                     method::Symbol = :spotsweeper) -> DataFrame

Identify likely technical artifacts (fiducials, tissue folds, tissue edges)
in `ds`. Returns a per-observation DataFrame with artifact probability scores.

Supported `method` values (planned):
- `:spotsweeper`  — SpotSweeper-style local outlier detection
- `:manual`       — returns empty scores, for interactive marking in SpatialViz
"""
function detect_artifacts(ds::SpatialDataset; method::Symbol = :spotsweeper)
    # TODO Phase 2: implement
    error("detect_artifacts: not yet implemented (Phase 2 deliverable)")
end

"""
    spatially_aware_qc(ds::SpatialDataset;
                       radius::Real = 100.0,
                       metrics::Vector{Symbol} = [:transcript_count, :gene_count])
        -> DataFrame

Compute spatially-aware QC metrics for each observation in `ds`, taking into
account the local neighbourhood within `radius` micrometres.
Returns a DataFrame with one row per observation and one column per metric.
"""
function spatially_aware_qc(
    ds::SpatialDataset;
    radius::Real = 100.0,
    metrics::Vector{Symbol} = [:transcript_count, :gene_count],
)
    # TODO Phase 2: implement
    error("spatially_aware_qc: not yet implemented (Phase 2 deliverable)")
end

end # module QualityControl
