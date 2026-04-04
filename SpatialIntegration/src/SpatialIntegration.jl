"""
    SpatialIntegration.jl

Cross-dataset harmonization and reference-mapping layer for the STX_DEV
spatial transcriptomics framework. Provides functions to align multiple
`SpatialDataset` objects into a common coordinate frame, correct
cross-experiment batch effects, and map query datasets onto annotated
reference atlases.
"""
module SpatialIntegration

using DataFrames
using SpatialOmicsBase
using SpatialAnalysis

# ---------------------------------------------------------------------------
# Exports — harmonization
# ---------------------------------------------------------------------------

export harmonize_datasets
export align_coordinate_systems
export standardize_features

# ---------------------------------------------------------------------------
# Exports — reference mapping
# ---------------------------------------------------------------------------

export integrate_reference
export transfer_labels
export deconvolute_spots

# ---------------------------------------------------------------------------
# Exports — multi-modal
# ---------------------------------------------------------------------------

export integrate_modalities
export joint_embedding

# ---------------------------------------------------------------------------
# Includes
# ---------------------------------------------------------------------------

include("harmonization/harmonize.jl")
include("harmonization/coordinate_alignment.jl")
include("harmonization/feature_standardization.jl")
include("reference/reference_mapping.jl")
include("reference/label_transfer.jl")
include("reference/deconvolution.jl")
include("multimodal/integration.jl")

end # module SpatialIntegration
