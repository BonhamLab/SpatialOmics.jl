"""
    SpatialStatistics

Spatial statistical methods: autocorrelation, hotspot detection, and
spatially-constrained clustering. Parallelised via ThreadsX.jl.
"""
module SpatialStatistics

using DataFrames
using Statistics
using ThreadsX
using SpatialOmicsBase

export spatial_autocorrelation, hotspot_detection, spatial_clustering

"""
    spatial_autocorrelation(ds::SpatialDataset,
                            feature::String;
                            method::Symbol = :morans_i,
                            n_permutations::Int = 999) -> NamedTuple

Compute a spatial autocorrelation statistic for `feature` across the spatial
observations in `ds`.

Returns a NamedTuple with fields `statistic`, `pvalue`, `zscore`.

Supported `method` values (planned):
- `:morans_i`   — Moran's I global autocorrelation
- `:gearys_c`   — Geary's C global autocorrelation
- `:getis_ord`  — Getis-Ord G* local statistic (hotspot variant)
"""
function spatial_autocorrelation(
    ds::SpatialDataset,
    feature::String;
    method::Symbol = :morans_i,
    n_permutations::Int = 999,
)
    # TODO Phase 2: implement; parallelise permutations with ThreadsX
    error("spatial_autocorrelation: not yet implemented (Phase 2 deliverable)")
end

"""
    hotspot_detection(ds::SpatialDataset,
                      feature::String;
                      fdr_threshold::Float64 = 0.05) -> DataFrame

Identify statistically significant spatial hotspots (high-expression regions)
for `feature`. Returns a per-observation DataFrame with columns:
`local_statistic`, `pvalue`, `fdr`, `is_hotspot`.
"""
function hotspot_detection(
    ds::SpatialDataset,
    feature::String;
    fdr_threshold::Float64 = 0.05,
)
    # TODO Phase 2: implement
    error("hotspot_detection: not yet implemented (Phase 2 deliverable)")
end

"""
    spatial_clustering(ds::SpatialDataset;
                       method::Symbol = :bayesspace,
                       n_clusters::Union{Int, Nothing} = nothing,
                       use_pca::Bool = true,
                       n_pcs::Int = 30) -> Vector{Int}

Cluster spatial observations into coherent tissue domains.
Returns a vector of integer cluster labels, one per observation.

Supported `method` values (planned):
- `:bayesspace`   — BayesSpace spatially-aware clustering (via PythonCall)
- `:leiden`       — Leiden graph clustering on spatial neighbourhood graph
- `:kmeans`       — k-means on PCA-reduced expression + spatial coordinates
"""
function spatial_clustering(
    ds::SpatialDataset;
    method::Symbol = :bayesspace,
    n_clusters::Union{Int, Nothing} = nothing,
    use_pca::Bool = true,
    n_pcs::Int = 30,
)
    # TODO Phase 2: implement
    error("spatial_clustering: not yet implemented (Phase 2 deliverable)")
end

end # module SpatialStatistics
