"""
    SpatialAnalysis.jl

Analysis core for the STX_DEV spatial transcriptomics framework.
Spatial statistics, QC, preprocessing, and cell interaction functions.
All public functions are pure (or explicitly documented as not pure),
enabling straightforward parallelisation via ThreadsX.jl and Distributed.jl.
"""
module SpatialAnalysis

using DataFrames
using Statistics
using ThreadsX
using SpatialOmicsBase

# ---------------------------------------------------------------------------
# Sub-modules (QC, Preprocessing, CellInteraction)
# ---------------------------------------------------------------------------

# TODO: All of these files should be direct in src/
# Subfolders unnecessary
include("qc/QualityControl.jl")
include("preprocessing/Preprocessing.jl")
include("interaction/CellInteraction.jl")

# TODO: Quit it with the namespaces. Packages don't need tons of submodules
using .QualityControl
using .Preprocessing
using .CellInteraction

# ---------------------------------------------------------------------------
# Statistics — included directly (no submodule namespace)
# ---------------------------------------------------------------------------

include("statistics/SpatialStatistics.jl")

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

# QualityControl
export filter_spots, detect_artifacts, spatially_aware_qc

# Preprocessing
export normalize_counts, correct_batch_effects, impute_missing, crop

# Statistics (owned by SpatialAnalysis directly)
export spatial_autocorrelation, hotspot_detection, spatial_clustering
export distances, classify_cells

# CellInteraction
export ligand_receptor_analysis, neighborhood_analysis, communication_score

# ---------------------------------------------------------------------------
# Top-level pipeline helpers
# ---------------------------------------------------------------------------

export run_qc_pipeline, default_qc_pipeline

include("pipeline.jl")

end # module SpatialAnalysis
