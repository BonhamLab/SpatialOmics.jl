# SpatialViz/src/recipes/shapes.jl
#
# Makie recipes for SpatialShapes:
#   poly!  — fills / outlines geometry
#   text!  — draws shape.name at the centroid of each geometry
#
# Geometries may be GeometryBasics.Polygon (cell boundaries, tissue outlines)
# or GeometryBasics.Circle (cell_circles, Visium spots). Both render via poly!.
# nothing entries (failed WKB decodes or HDF5 shapes) are silently skipped.
#
# Usage:
#   poly!(ax, shp; color=:transparent, strokecolor=:white, strokewidth=0.5)
#   text!(ax, shp; fontsize=8, align=(:center,:center))
#   text!(ax, shp["roi_1"])

import SpatialOmicsBase: _geom_centroid

const _SpatialGeom = Union{GeometryBasics.Polygon, GeometryBasics.Circle}

# ---------------------------------------------------------------------------
# poly! — geometry fill / outline
# ---------------------------------------------------------------------------

function Makie.convert_arguments(P::Type{<:Poly}, v::SpatialElementView{<:SpatialShapes})
    return Makie.convert_arguments(P, collect(v))
end

"""
    Makie.convert_arguments(P, shp::SpatialShapes)

Convert a `SpatialShapes` to poly arguments. Accepts `Polygon` and `Circle`
geometries; `nothing` entries are silently dropped.
"""
function Makie.convert_arguments(P::Type{<:Poly}, shp::SpatialShapes)
    geoms = _SpatialGeom[s.geometry for s in shp.shapes if s.geometry isa _SpatialGeom]
    return Makie.convert_arguments(P, geoms)
end

# ---------------------------------------------------------------------------
# text! — shape names at centroids
# ---------------------------------------------------------------------------

function _shape_label_args(shp::SpatialShapes)
    pts    = Point2f[]
    labels = String[]
    for s in shp.shapes
        cxy = _geom_centroid(s.geometry)
        cxy === nothing && continue
        push!(pts,    Point2f(cxy[1], cxy[2]))
        push!(labels, s.name)
    end
    return pts, labels
end

"""
    text!(ax, shp::SpatialShapes; kwargs...)

Draw each shape's `name` at the centroid of its geometry. Shapes with no
computable centroid (e.g. `nothing` geometry) are skipped. All Makie `text!`
keyword arguments (`fontsize`, `align`, `color`, …) are forwarded.
"""
function Makie.text!(ax, shp::SpatialShapes; kwargs...)
    pts, labels = _shape_label_args(shp)
    isempty(pts) && return nothing
    return Makie.text!(ax, pts; text = labels, kwargs...)
end

"""
    text!(ax, s::SpatialShape; kwargs...)

Draw `s.name` at the centroid of `s.geometry`.
"""
function Makie.text!(ax, s::SpatialShape; kwargs...)
    cxy = _geom_centroid(s.geometry)
    cxy === nothing && return nothing
    return Makie.text!(ax, [Point2f(cxy[1], cxy[2])]; text = [s.name], kwargs...)
end

function Makie.text!(ax, v::SpatialElementView{<:SpatialShapes}; kwargs...)
    return Makie.text!(ax, collect(v); kwargs...)
end
