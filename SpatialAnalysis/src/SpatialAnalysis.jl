"""
    SpatialAnalysis.jl

Analysis core for the STX_DEV spatial transcriptomics framework.
Provides four functional sub-modules: QualityControl, Preprocessing,
SpatialStatistics, and CellInteraction. All public functions are pure (or
explicitly documented as not pure), enabling straightforward parallelisation
via ThreadsX.jl and Distributed.jl.
"""
module SpatialAnalysis

using DataFrames
using Statistics
using ThreadsX
using SpatialOmicsBase

# ---------------------------------------------------------------------------
# Sub-modules
# ---------------------------------------------------------------------------

include("qc/QualityControl.jl")
include("preprocessing/Preprocessing.jl")
include("statistics/SpatialStatistics.jl")
include("interaction/CellInteraction.jl")

# Bring sub-module exports into the SpatialAnalysis namespace so callers can
# do either `using SpatialAnalysis` or `using SpatialAnalysis.QualityControl`.
using .QualityControl
using .Preprocessing
using .SpatialStatistics
using .CellInteraction

# ---------------------------------------------------------------------------
# Re-export everything from sub-modules
# ---------------------------------------------------------------------------

# QualityControl
export filter_spots, detect_artifacts, spatially_aware_qc

# Preprocessing
export normalize_counts, correct_batch_effects, impute_missing

# SpatialStatistics
export spatial_autocorrelation, hotspot_detection, spatial_clustering

# CellInteraction
export ligand_receptor_analysis, neighborhood_analysis, communication_score

# ---------------------------------------------------------------------------
# Top-level pipeline helpers
# ---------------------------------------------------------------------------

export run_qc_pipeline, default_qc_pipeline

include("pipeline.jl")

end # module SpatialAnalysis
