# SpatialOmics/src/types/coordinate_systems.jl
# Coordinate system definitions following the SpatialData spec.

"""
    CoordinateSystem

Describes a named coordinate space with physical units.

Fields
------
- `name`  : String — e.g. "global", "tissue_section_1"
- `axes`  : Vector{String} — ordered axis names, e.g. ["x", "y", "z"]
- `units` : Vector{String} — corresponding units, e.g. ["µm", "µm"]
"""
struct CoordinateSystem
    name::String
    axes::Vector{String}
    units::Vector{String}
end

# Convenience constructor for 2-D physical space (micrometres)
CoordinateSystem(name::String) = CoordinateSystem(name, ["x", "y"], ["µm", "µm"])

function Base.show(io::IO, cs::CoordinateSystem)
    axes_str  = "[" * join(cs.axes, ", ") * "]"
    units_str = "[" * join(repr.(cs.units), ", ") * "]"
    print(io, "CoordinateSystem($(repr(cs.name)), axes=$(axes_str), units=$(units_str))")
end
