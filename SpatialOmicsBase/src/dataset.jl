# SpatialOmics/src/types/dataset.jl
# The top-level SpatialDataset container.

"""
    SpatialDataset{T<:AbstractFloat}

The primary data container for the STX_DEV framework, following the SpatialData
specification. All spatial elements are stored in named dictionaries; coordinate
systems and transformations are tracked explicitly so multi-modal alignment is
a first-class operation.

Type parameter `T` is the floating-point precision used for coordinates and
image pixel values (typically Float32 or Float64).

Fields
------
- `images`              : Dict{String, SpatialImage{T}}
- `labels`              : Dict{String, SpatialLabels}
- `points`              : Dict{String, SpatialPoints{T}}
- `shapes`              : Dict{String, SpatialShapes}
- `tables`              : Dict{String, SpatialTable}
- `coordinate_systems`  : Dict{String, CoordinateSystem}
- `transformations`     : Dict{String, Transformation}
- `metadata`            : Dict{String, Any}
"""
mutable struct SpatialDataset{T<:AbstractFloat}
    images::Dict{String, SpatialImage}       # SpatialImage{T} for any T (uint8, float32, …)
    labels::Dict{String, SpatialLabels}
    points::Dict{String, SpatialPoints{T}}
    shapes::Dict{String, SpatialShapes}
    tables::Dict{String, SpatialTable}
    coordinate_systems::Dict{String, CoordinateSystem}
    transformations::Dict{String, Transformation}
    metadata::Dict{String, Any}
end

"""
    spatial_dataset(; T=Float32, metadata=Dict{String,Any}()) -> SpatialDataset{T}

Construct an empty `SpatialDataset` with the given floating-point precision.
"""
function spatial_dataset(;
    T::Type{<:AbstractFloat} = Float32,
    metadata::Dict{String,Any} = Dict{String,Any}(),
)
    return SpatialDataset{T}(
        Dict{String, SpatialImage}(),
        Dict{String, SpatialLabels}(),
        Dict{String, SpatialPoints{T}}(),
        Dict{String, SpatialShapes}(),
        Dict{String, SpatialTable}(),
        Dict{String, CoordinateSystem}(),
        Dict{String, Transformation}(),
        metadata,
    )
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------

function Base.show(io::IO, ds::SpatialDataset{T}) where T
    println(io, "SpatialDataset{$T}")
    _show_element_section(io, "images",  ds.images,  _summary_image)
    _show_element_section(io, "labels",  ds.labels,  _summary_labels)
    _show_element_section(io, "points",  ds.points,  _summary_points)
    _show_element_section(io, "shapes",  ds.shapes,  _summary_shapes)
    _show_element_section(io, "tables",  ds.tables,  _summary_table)
    if !isempty(ds.coordinate_systems)
        cs_names = join(keys(ds.coordinate_systems), ", ")
        println(io, "  coord_systems ($(length(ds.coordinate_systems))): $cs_names")
    end
end

function _show_element_section(io, label, dict, summarise)
    isempty(dict) && return
    n = length(dict)
    pad = " " ^ (8 - length(label))   # align colons
    print(io, "  $label$pad($n):")
    for (name, el) in dict
        print(io, "  $name [$(summarise(el))]")
    end
    println(io)
end

_summary_image(img::SpatialImage)    = join(size(img.data), "×") * " $(eltype(img.data))"
_summary_labels(lbl::SpatialLabels)  = join(size(lbl.data), "×") * " $(eltype(lbl.data))"
_summary_points(pts::SpatialPoints)  = "$(size(pts.coordinates, 1)) pts × $(ncol(pts.features)) cols"
_summary_shapes(shp::SpatialShapes)  = "$(length(shp.geometries)) shapes × $(ncol(shp.features)) cols"
_summary_table(tbl::SpatialTable)    = "$(size(tbl.data, 1)) obs × $(size(tbl.data, 2)) vars"

# ---------------------------------------------------------------------------
# Indexing — ds[name] = element  /  ds[name, ElementType]
# ---------------------------------------------------------------------------

Base.setindex!(ds::SpatialDataset, img::SpatialImage,   k::String) = (ds.images[k]  = img;  ds)
Base.setindex!(ds::SpatialDataset, pts::SpatialPoints,  k::String) = (ds.points[k]  = pts;  ds)
Base.setindex!(ds::SpatialDataset, lbl::SpatialLabels,  k::String) = (ds.labels[k]  = lbl;  ds)
Base.setindex!(ds::SpatialDataset, shp::SpatialShapes,  k::String) = (ds.shapes[k]  = shp;  ds)
Base.setindex!(ds::SpatialDataset, tbl::SpatialTable,   k::String) = (ds.tables[k]  = tbl;  ds)

"""
    images(ds)              -> Dict{String, SpatialImage}
    images(ds, name)        -> SpatialImage
    images(roi)             -> ViewDict (scoped to roi.extent)
    images(roi, name)       -> SpatialElementView{SpatialImage}
"""
images(ds::SpatialDataset)              = ds.images
images(ds::SpatialDataset, k::String)   = ds.images[k]

"""
    labels(ds)              -> Dict{String, SpatialLabels}
    labels(ds, name)        -> SpatialLabels
"""
labels(ds::SpatialDataset)              = ds.labels
labels(ds::SpatialDataset, k::String)   = ds.labels[k]

"""
    points(ds)              -> Dict{String, SpatialPoints}
    points(ds, name)        -> SpatialPoints
"""
points(ds::SpatialDataset)              = ds.points
points(ds::SpatialDataset, k::String)   = ds.points[k]

"""
    shapes(ds)              -> Dict{String, SpatialShapes}
    shapes(ds, name)        -> SpatialShapes
"""
shapes(ds::SpatialDataset)              = ds.shapes
shapes(ds::SpatialDataset, k::String)   = ds.shapes[k]

"""
    tables(ds)              -> Dict{String, SpatialTable}
    tables(ds, name)        -> SpatialTable
"""
tables(ds::SpatialDataset)              = ds.tables
tables(ds::SpatialDataset, k::String)   = ds.tables[k]

"""
    metadata(ds)            -> Dict{String, Any}
"""
metadata(ds::SpatialDataset)            = ds.metadata
