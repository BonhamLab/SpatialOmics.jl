# SpatialIO/src/formats/geojson.jl
# Export SpatialShapes as GeoJSON for use in GIS tools and Napari.

"""
    write(shapes::SpatialShapes, path::String, ::GeoJSON)

Write `shapes` to a GeoJSON FeatureCollection file at `path`.
Each geometry is serialised as a GeoJSON Feature with properties drawn
from `shapes.features`.
"""
function write(shapes::SpatialShapes, path::String, ::GeoJSON)
    # TODO Phase 1: implement with JSON.jl
    error("write GeoJSON: not yet implemented (Phase 1 deliverable)")
end
