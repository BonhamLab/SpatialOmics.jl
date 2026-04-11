# SpatialOmics/src/transforms/affine.jl
# Implementations of transform_coordinates, compose_transformations, and
# invert_transformation for AffineTransformation.
# TODO: Shouldn't these use JuliaImages ecosystem?

function transform_coordinates(t::IdentityTransformation, coords::Matrix{<:AbstractFloat})
    return copy(coords)
end

function transform_coordinates(t::AffineTransformation, coords::Matrix{<:AbstractFloat})
    # TODO Phase 1: apply t.matrix to homogeneous coords
    error("AffineTransformation.transform_coordinates: not yet implemented")
end

function compose_transformations(a::Transformation, b::Transformation)
    # TODO Phase 1: return a new transformation equivalent to applying a then b
    error("compose_transformations: not yet implemented")
end

function invert_transformation(t::IdentityTransformation)
    return IdentityTransformation()
end

function invert_transformation(t::AffineTransformation)
    # TODO Phase 1: invert augmented affine matrix
    error("AffineTransformation.invert_transformation: not yet implemented")
end
