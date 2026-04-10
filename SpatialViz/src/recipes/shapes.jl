# SpatialViz/src/recipes/shapes.jl
#
# convert_arguments for SpatialShapes → poly
#
# Geometries may be GeometryBasics.Polygon (cell boundaries, tissue outlines)
# or GeometryBasics.Circle (cell_circles, Visium spots). Both render via poly!.
# nothing entries (failed WKB decodes) are silently dropped.
#
# Usage:
#   poly(shp)
#   poly!(ax, shp; color=:transparent, strokecolor=:white, strokewidth=0.5)

const _SpatialGeom = Union{GeometryBasics.Polygon, GeometryBasics.Circle}

# View dispatch — crop is applied lazily at render time.
function Makie.convert_arguments(P::Type{<:Poly}, v::SpatialElementView{<:SpatialShapes})
    return Makie.convert_arguments(P, crop(v.parent, v.extent))
end

"""
    Makie.convert_arguments(P, shp::SpatialShapes)

Convert a `SpatialShapes` to poly arguments. Accepts `Polygon` and `Circle`
geometries; `nothing` entries (failed WKB decodes) are silently dropped.
"""
function Makie.convert_arguments(P::Type{<:Poly}, shp::SpatialShapes)
    geoms = _SpatialGeom[g for g in shp.geometries if g isa _SpatialGeom]
    return Makie.convert_arguments(P, geoms)
end
