"""
    SpatialOmicsBase.jl

Lightweight core types package for the STX_DEV spatial transcriptomics framework.
Defines SpatialDataset and all spatial element types following the SpatialData
specification. Contains no I/O dependencies — all file reading/writing lives in
SpatialIO.jl, which depends on this package.

Other packages (SpatialIO, SpatialAnalysis, SpatialViz, SpatialIntegration) depend
on SpatialOmicsBase. End users depend on the SpatialOmics.jl umbrella.
"""
module SpatialOmicsBase

using DiskArrays
using DataFrames

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

export SpatialDataset
export SpatialElement
export SpatialImage, SpatialPoints, SpatialLabels, SpatialShapes, SpatialTable
export CoordinateSystem, Transformation
export AffineTransformation, IdentityTransformation

# Construction helpers
export spatial_dataset

export images, labels, points, shapes, tables, metadata
export channels, channels!

# Coordinate utilities
export transform_coordinates, compose_transformations, invert_transformation

# Metadata helpers
export get_metadata, set_metadata!

# Spatial query primitives
export SpatialExtent
export SpatialElementView, SpatialDatasetView
export extent, intersects, crop

# Pyramid sampler (for Makie.Resampler integration)
export ImagePyramidSampler

# ---------------------------------------------------------------------------
# Includes (one file per logical group)
# ---------------------------------------------------------------------------

include("types/elements.jl")
include("types/coordinate_systems.jl")
include("types/transformations.jl")
include("types/dataset.jl")
include("types/view.jl")
include("types/pyramid.jl")
include("transforms/affine.jl")
include("utils/metadata.jl")
include("utils/chunking.jl")

end # module SpatialOmicsBase
