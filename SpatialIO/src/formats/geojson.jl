# SpatialIO/src/formats/geojson.jl
# Export SpatialShapes as GeoJSON for use in GIS tools and Napari.

"""
    write_geojson(shapes::SpatialShapes, path::String)

Write `shapes` to a GeoJSON FeatureCollection file at `path`.
Each geometry is serialised as a GeoJSON Feature with properties drawn
from `shapes.features`.
"""
function write_geojson(shapes::SpatialShapes, path::String)
    # TODO Phase 1: implement with JSON3.jl
    error("write_geojson: not yet implemented (Phase 1 deliverable)")
end
