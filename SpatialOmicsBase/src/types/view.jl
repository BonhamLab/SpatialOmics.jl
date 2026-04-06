# SpatialOmicsBase/src/types/view.jl
# Spatial query primitives: bounding boxes and scoped views of a SpatialDataset.
#
# These types are pure coordinate / data structures — no rendering dependency.
# They underpin the three visualization scopes:
#   Slide  → SpatialView over the global coordinate system, no extent filter
#   FOV    → SpatialView in a local coordinate system, no extent filter
#   Region → SpatialView with an SpatialExtent to crop to a subregion

# ---------------------------------------------------------------------------
# SpatialExtent
# ---------------------------------------------------------------------------

"""
    SpatialExtent{T<:AbstractFloat}

An axis-aligned bounding box in a named coordinate system.

# Fields
- `xmin`, `xmax` : x (column) bounds
- `ymin`, `ymax` : y (row) bounds
- `coordinate_system` : name of the coordinate system these bounds live in

# Constructors
```julia
SpatialExtent(0, 100, 0, 200)                         # Float64, cs="global"
SpatialExtent(0, 100, 0, 200, "slide_1")              # named coordinate system
SpatialExtent{Float32}(0, 100, 0, 200, "global")      # explicit element type
SpatialExtent((0.0, 100.0), (0.0, 200.0), "slide_1")  # tuple form
```
"""
struct SpatialExtent{T<:AbstractFloat}
    xmin::T
    xmax::T
    ymin::T
    ymax::T
    coordinate_system::String
end

# Convenience constructors — accept plain Real without explicit type suffix.
# Untyped form defaults to Float64. Typed form SpatialExtent{T}(x,y,z,w) works
# via the struct's inner constructor (which uses convert internally) — no outer
# parametric wrapper needed (it would shadow the inner constructor and recurse).
SpatialExtent(xmin::Real, xmax::Real, ymin::Real, ymax::Real, cs::String="global") =
    SpatialExtent{Float64}(Float64(xmin), Float64(xmax), Float64(ymin), Float64(ymax), cs)

# Tuple-based convenience constructors
SpatialExtent(x::Tuple{Real,Real}, y::Tuple{Real,Real}, cs::String="global") =
    SpatialExtent{Float64}(Float64(x[1]), Float64(x[2]), Float64(y[1]), Float64(y[2]), cs)

SpatialExtent{T}(x::Tuple{Real,Real}, y::Tuple{Real,Real}, cs::String="global") where {T<:AbstractFloat} =
    SpatialExtent{T}(T(x[1]), T(x[2]), T(y[1]), T(y[2]), cs)

function Base.show(io::IO, ext::SpatialExtent{T}) where T
    print(io, "SpatialExtent{$T}(",
          "x=[$(ext.xmin), $(ext.xmax)], ",
          "y=[$(ext.ymin), $(ext.ymax)], ",
          "cs=\"$(ext.coordinate_system)\")")
end

"""Width of the extent in x."""
width(ext::SpatialExtent)  = ext.xmax - ext.xmin

"""Height of the extent in y."""
height(ext::SpatialExtent) = ext.ymax - ext.ymin

"""Test whether point `(x, y)` lies within the extent (inclusive)."""
function Base.in(pt::Tuple{Real,Real}, ext::SpatialExtent)
    return ext.xmin <= pt[1] <= ext.xmax && ext.ymin <= pt[2] <= ext.ymax
end

"""
    intersects(a::SpatialExtent, b::SpatialExtent) -> Bool

Return `true` if `a` and `b` overlap. Extents in different coordinate systems
never intersect.
"""
function intersects(a::SpatialExtent, b::SpatialExtent)
    a.coordinate_system == b.coordinate_system || return false
    return a.xmin <= b.xmax && a.xmax >= b.xmin &&
           a.ymin <= b.ymax && a.ymax >= b.ymin
end

# ---------------------------------------------------------------------------
# extent() — compute the bounding box of a spatial element
# ---------------------------------------------------------------------------

"""
    extent(pts::SpatialPoints) -> SpatialExtent

Return the axis-aligned bounding box of the point coordinates.
"""
function extent(pts::SpatialPoints{T})::SpatialExtent{T} where T
    cs = get(pts.metadata, "coordinate_system", "global")
    return SpatialExtent{T}(
        minimum(pts.coordinates[:, 1]), maximum(pts.coordinates[:, 1]),
        minimum(pts.coordinates[:, 2]), maximum(pts.coordinates[:, 2]),
        cs,
    )
end

"""
    extent(img::SpatialImage) -> SpatialExtent

Return the pixel-coordinate bounding box of the image (1-based, last two dims
are y × x following the OME-Zarr c,y,x convention).
"""
function extent(img::SpatialImage{T})::SpatialExtent{T} where T
    sz = size(img.data)
    ny, nx = sz[end-1], sz[end]
    cs = get(img.metadata, "coordinate_system", "global")
    return SpatialExtent{T}(T(1), T(nx), T(1), T(ny), cs)
end

"""
    extent(lbl::SpatialLabels) -> SpatialExtent{Float32}

Return the pixel bounding box of the label image.
"""
function extent(lbl::SpatialLabels)::SpatialExtent{Float32}
    sz = size(lbl.data)
    ny, nx = sz[end-1], sz[end]
    cs = get(lbl.metadata, "coordinate_system", "global")
    return SpatialExtent{Float32}(1f0, Float32(nx), 1f0, Float32(ny), cs)
end

# ---------------------------------------------------------------------------
# crop() — filter elements to an extent
# ---------------------------------------------------------------------------

"""
    crop(pts::SpatialPoints, ext::SpatialExtent) -> SpatialPoints

Return a new `SpatialPoints` containing only the points within `ext`.
"""
function crop(pts::SpatialPoints{T}, ext::SpatialExtent) where T
    x = pts.coordinates[:, 1]
    y = pts.coordinates[:, 2]
    mask = (x .>= ext.xmin) .& (x .<= ext.xmax) .&
           (y .>= ext.ymin) .& (y .<= ext.ymax)
    return SpatialPoints(pts.coordinates[mask, :], pts.features[mask, :], pts.metadata)
end

"""
    crop(shp::SpatialShapes, ext::SpatialExtent) -> SpatialShapes

Return shapes whose bounding boxes overlap `ext`. (Fast pre-filter using
centroid or first vertex; exact intersection is left to the caller.)
"""
function crop(shp::SpatialShapes, ext::SpatialExtent)
    # Keep all shapes for now — full geometric intersection deferred to Phase 2
    return shp
end

# ---------------------------------------------------------------------------
# SpatialView
# ---------------------------------------------------------------------------

"""
    SpatialView

A scoped query into a `SpatialDataset`: a target coordinate system, an optional
spatial extent for cropping, and an optional allow-list of element keys.

# Constructors
```julia
SpatialView("global")                          # entire slide in global CS
SpatialView("global", ext)                     # crop to SpatialExtent
SpatialView("global", ext, ["morphology", "transcripts"])  # select layers
```
"""
struct SpatialView
    coordinate_system::String
    extent::Union{SpatialExtent, Nothing}
    element_keys::Union{Vector{String}, Nothing}
end

SpatialView(cs::String) = SpatialView(cs, nothing, nothing)
SpatialView(cs::String, ext::SpatialExtent) = SpatialView(cs, ext, nothing)

function Base.show(io::IO, v::SpatialView)
    print(io, "SpatialView(cs=\"$(v.coordinate_system)\"")
    v.extent !== nothing && print(io, ", extent=$(v.extent)")
    v.element_keys !== nothing && print(io, ", keys=$(v.element_keys)")
    print(io, ")")
end
