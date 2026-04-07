# SpatialViz/src/recipes/points.jl
#
# convert_arguments for SpatialPoints → scatter
#
# Usage:
#   scatter(pts)                         # plain scatter, no color
#   scatter(pts; color=pts.features.gene)  # color by a feature column

# View dispatch — crop is applied lazily at render time.
function Makie.convert_arguments(P::Type{<:Scatter}, v::SpatialElementView{<:SpatialPoints})
    return Makie.convert_arguments(P, crop(v.parent, v.extent))
end

function Makie.convert_arguments(P::Type{<:Lines}, v::SpatialElementView{<:SpatialPoints})
    return Makie.convert_arguments(P, crop(v.parent, v.extent))
end

"""
    Makie.convert_arguments(P, pts::SpatialPoints)

Convert a `SpatialPoints` to scatter arguments: a vector of `Point2f`
constructed from the first two coordinate columns (x = col 1, y = col 2).
"""
function Makie.convert_arguments(P::Type{<:Scatter}, pts::SpatialPoints)
    coords = pts.coordinates
    points = [Point2f(coords[i, 1], coords[i, 2]) for i in axes(coords, 1)]
    return Makie.convert_arguments(P, points)
end

# Also register for Lines so users can draw point trajectories if needed.
function Makie.convert_arguments(P::Type{<:Lines}, pts::SpatialPoints)
    coords = pts.coordinates
    points = [Point2f(coords[i, 1], coords[i, 2]) for i in axes(coords, 1)]
    return Makie.convert_arguments(P, points)
end
