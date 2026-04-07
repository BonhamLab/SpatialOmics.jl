# SpatialViz/src/recipes/shapes.jl
#
# convert_arguments for SpatialShapes → poly
#
# Geometries are expected to be GeometryBasics.Polygon values (as decoded by
# SpatialIO's WKB reader). Non-Polygon entries (e.g. nothing for failed WKB
# decodes) are silently dropped.
#
# Usage:
#   poly(shp)
#   poly!(ax, shp; color=:transparent, strokecolor=:white, strokewidth=0.5)

# View dispatch — crop is applied lazily at render time.
function Makie.convert_arguments(P::Type{<:Poly}, v::SpatialElementView{<:SpatialShapes})
    return Makie.convert_arguments(P, crop(v.parent, v.extent))
end

"""
    Makie.convert_arguments(P, shp::SpatialShapes)

Convert a `SpatialShapes` to poly arguments: a `Vector{GeometryBasics.Polygon}`
containing only the entries that decoded successfully (non-nothing, Polygon type).
"""
function Makie.convert_arguments(P::Type{<:Poly}, shp::SpatialShapes)
    polys = GeometryBasics.Polygon[g for g in shp.geometries
                                   if g isa GeometryBasics.Polygon]
    return Makie.convert_arguments(P, polys)
end
