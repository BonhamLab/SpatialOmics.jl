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
using GeometryBasics
using GeometryOps

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
export region, region_key, instance_key

# Coordinate utilities
export transform_coordinates, compose_transformations, invert_transformation

# Metadata helpers
export get_metadata, set_metadata!

# Spatial query primitives
export SpatialExtent
export SpatialElementView, SpatialDatasetView
export extent, intersects, crop, subset

# Pyramid sampler (for Makie.Resampler integration)
export ImagePyramidSampler

# ---------------------------------------------------------------------------
# Includes (one file per logical group)
# ---------------------------------------------------------------------------

include("elements.jl")
include("coordinate_systems.jl")
include("transformations.jl")
include("dataset.jl")
include("view.jl")
include("pyramid.jl")
include("affine.jl")
include("metadata.jl")
include("chunking.jl")

end # module SpatialOmicsBase
