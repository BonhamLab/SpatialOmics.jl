# SpatialOmicsBase/src/affine.jl
# Coordinate transform implementations using CoordinateTransformations.jl.
#
# AffineTransformation wraps an (D+1)×(D+1) augmented matrix.
# We delegate to CoordinateTransformations.AffineMap for the actual maths so
# that downstream code can compose and invert without reimplementing linear algebra.

using CoordinateTransformations
using LinearAlgebra

"""
    transform_coordinates(t::IdentityTransformation, coords::AbstractMatrix{<:AbstractFloat})
        -> AbstractMatrix

Return a copy of `coords` (identity transform is a no-op).
"""
function transform_coordinates(t::IdentityTransformation, coords::AbstractMatrix{<:AbstractFloat})
    return copy(coords)
end

"""
    transform_coordinates(t::AffineTransformation, coords::AbstractMatrix{<:AbstractFloat})
        -> Matrix{Float64}

Apply the affine transformation `t` to `coords` (N×D matrix, rows = points).
Returns an N×D `Matrix{Float64}`.
"""
function transform_coordinates(t::AffineTransformation, coords::AbstractMatrix{<:AbstractFloat})
    D = size(t.matrix, 1) - 1   # spatial dimension
    A = t.matrix[1:D, 1:D]      # linear part
    b = t.matrix[1:D, end]      # translation
    amap = CoordinateTransformations.AffineMap(A, b)
    # Apply row-wise: each row of coords is a D-dimensional point
    return Matrix{Float64}(reduce(vcat, [amap(Float64.(coords[i, :]))' for i in axes(coords, 1)]))
end

"""
    compose_transformations(a::Transformation, b::Transformation) -> Transformation

Return a new transformation equivalent to applying `a` then `b`.
`IdentityTransformation` is the neutral element on both sides.
"""
function compose_transformations(a::IdentityTransformation, b::Transformation)
    return b
end

function compose_transformations(a::Transformation, b::IdentityTransformation)
    return a
end

function compose_transformations(a::AffineTransformation, b::AffineTransformation)
    M = b.matrix * a.matrix   # augmented matrix composition: apply a first, then b
    return AffineTransformation(M, a.input_coordinate_system, b.output_coordinate_system)
end

"""
    invert_transformation(t::IdentityTransformation) -> IdentityTransformation

Return the identity (self-inverse).
"""
function invert_transformation(t::IdentityTransformation)
    return IdentityTransformation()
end

"""
    invert_transformation(t::AffineTransformation) -> AffineTransformation

Return a new `AffineTransformation` whose augmented matrix is the inverse of `t.matrix`.
Throws `LinearAlgebra.SingularException` if the matrix is not invertible.
"""
function invert_transformation(t::AffineTransformation)
    M_inv = inv(t.matrix)
    return AffineTransformation(M_inv, t.output_coordinate_system, t.input_coordinate_system)
end
