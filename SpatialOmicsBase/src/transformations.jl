# SpatialOmics/src/types/transformations.jl
# Spatial transformations between coordinate systems.

"""
    Transformation

Abstract supertype for all coordinate transformations.
"""
abstract type Transformation end

"""
    IdentityTransformation

The do-nothing transformation. Useful as a default when elements already live
in the target coordinate system.
"""
struct IdentityTransformation <: Transformation end

"""
    AffineTransformation

A general affine (linear + translation) transformation represented by an
augmented matrix M of size (D+1)×(D+1).

Fields
------
- `matrix`          : Matrix{Float64} — the (D+1)×(D+1) augmented affine matrix
- `input_coordinate_system`  : String
- `output_coordinate_system` : String
"""
struct AffineTransformation <: Transformation
    matrix::Matrix{Float64}
    input_coordinate_system::String
    output_coordinate_system::String
end

# Stub — implementation lives in transforms/affine.jl
function transform_coordinates end
function compose_transformations end
function invert_transformation end
