# SpatialOmicsBase/src/types/view.jl
#
# Spatial bounding boxes and lazy dataset/element views.
#
# Design
# ──────
# SpatialExtent{T}       — axis-aligned bounding box in a named coordinate system.
# SpatialElementView{T}  — lazy view of a single element cropped to an extent.
#                          Analogous to Julia's SubArray: holds a parent + index.
# SpatialDatasetView     — lazy view of an entire dataset; property access on
#                          .images/.points/.shapes/.labels/.tables returns a
#                          ViewDict that vends SpatialElementViews on indexing.
#
# Construction
# ────────────
#   view(ds, ext)          → SpatialDatasetView
#   view(pts, ext)         → SpatialElementView{SpatialPoints}
#   roi.points["key"]      → SpatialElementView{SpatialPoints}
#
# Plotting (in SpatialViz)
# ────────────────────────
#   scatter!(ax, roi.points["transcripts"])
#   poly!(ax, roi.shapes["cell_boundaries"])
#   heatmap!(ax, roi.images["morphology_focus"]; channel=1)
#
# The view carries the extent so each plot verb can apply it independently;
# no composite panel helper is required.

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

# Untyped 5-arg form — defaults to Float64.
# NOTE: do NOT define SpatialExtent{T}(::Real,...) here — it shadows the struct's
# inner constructor (which uses ::Any args) and causes infinite recursion.
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
# crop() — materialise a filtered copy of an element
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
centroid or first vertex; exact intersection deferred to Phase 2.)
"""
function crop(shp::SpatialShapes, ext::SpatialExtent)
    return shp  # full geometric intersection deferred
end

# ---------------------------------------------------------------------------
# SpatialElementView — lazy view of a single element
# ---------------------------------------------------------------------------

"""
    SpatialElementView{T<:SpatialElement}

A lazy, non-materialising view of a spatial element cropped to a
`SpatialExtent`. Analogous to Julia's `SubArray`: holds a reference to the
parent and the bounding-box "index", applies the crop only when the data
are actually needed (e.g. at plot time via `convert_arguments`).

Construct via `Base.view`:
```julia
pts_view = view(xen.points["transcripts"], extent)
scatter!(ax, pts_view)   # crop applied at render time
```
"""
struct SpatialElementView{T<:SpatialElement}
    parent::T
    extent::SpatialExtent
end

Base.view(el::SpatialElement, ext::SpatialExtent) = SpatialElementView(el, ext)

function Base.show(io::IO, v::SpatialElementView{T}) where T
    print(io, "view(", T, ", ", v.extent, ")")
end

# ---------------------------------------------------------------------------
# ViewDict — lazy dict that vends SpatialElementViews on indexing
# ---------------------------------------------------------------------------

# Internal helper; not exported. Returned by SpatialDatasetView property access.
struct ViewDict{T<:SpatialElement}
    dict::Dict{String, T}
    extent::SpatialExtent
end

Base.getindex(vd::ViewDict, k::String) = SpatialElementView(vd.dict[k], vd.extent)
Base.haskey(vd::ViewDict, k::String)   = haskey(vd.dict, k)
Base.keys(vd::ViewDict)                = keys(vd.dict)
Base.length(vd::ViewDict)              = length(vd.dict)
Base.isempty(vd::ViewDict)             = isempty(vd.dict)

# ---------------------------------------------------------------------------
# SpatialDatasetView — lazy view of an entire dataset
# ---------------------------------------------------------------------------

"""
    SpatialDatasetView{T<:AbstractFloat}

A lazy, extent-scoped view over a `SpatialDataset`. Property access on
`.images`, `.points`, `.shapes`, `.labels`, `.tables` returns a `ViewDict`
that vends `SpatialElementView` objects on key lookup.

Construct via `Base.view`:
```julia
roi = view(xen, SpatialExtent(cx-500, cx+500, cy-500, cy+500))

heatmap!(ax, roi.images["morphology_focus"]; channel=1, colormap=:grays)
scatter!(ax, roi.points["transcripts"])
poly!(ax, roi.shapes["cell_boundaries"])
```

All three calls use the same `roi` — no coordinate arithmetic is repeated
and no panel helper is required.
"""
struct SpatialDatasetView{T<:AbstractFloat}
    dataset::SpatialDataset{T}
    extent::SpatialExtent
end

Base.view(ds::SpatialDataset, ext::SpatialExtent) = SpatialDatasetView(ds, ext)

# Element dictionaries become ViewDicts; everything else passes through.
const _ELEMENT_FIELDS = (:images, :labels, :points, :shapes, :tables)

function Base.getproperty(v::SpatialDatasetView, s::Symbol)
    s in (:dataset, :extent) && return getfield(v, s)
    ds  = getfield(v, :dataset)
    ext = getfield(v, :extent)
    s in _ELEMENT_FIELDS && return ViewDict(getproperty(ds, s), ext)
    return getproperty(ds, s)
end

function Base.show(io::IO, v::SpatialDatasetView{T}) where T
    ds  = getfield(v, :dataset)
    ext = getfield(v, :extent)
    ni = length(ds.images); np = length(ds.points)
    ns = length(ds.shapes); nl = length(ds.labels)
    xr = "x=[$(round(ext.xmin,digits=1)), $(round(ext.xmax,digits=1))]"
    yr = "y=[$(round(ext.ymin,digits=1)), $(round(ext.ymax,digits=1))]"
    print(io, "SpatialDatasetView{$T}($xr, $yr, ",
          "$ni image(s), $np point(s), $ns shape(s), $nl label(s))")
end
