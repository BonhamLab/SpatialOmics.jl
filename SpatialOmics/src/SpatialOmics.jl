module SpatialOmics

using GeometryBasics: Point, Point2f, Polygon, AbstractGeometry
using GeoInterface
using Tables
using StaticArrays
using OrderedCollections
using GeometryOps
using Zarr
using JSON3
using ImageBase: restrict

export
    # Coordinate systems
    CoordinateSystem,
    # Transformations
    AbstractTransformation,
    Identity, Affine, Sequence,
    translation, scaling, rotation, flip_y, compose,
    apply, apply!, resolve_transform,
    # Dataset
    BackingStore, SpatialDataset,
    coord_systems, transform,
    with_dataset, keep!,
    # Elements
    SpatialPoints, SpatialShapes,
    points, shapes,
    coords, features, feature_ids, coord_system,
    geometries, bbox, instance_ids,
    # Views
    SpatialExtent, SpatialROI,
    SpatialElementView, SpatialDatasetView,
    geometry,
    # Images + Labels
    SpatialImage, SpatialLabels,
    data, nchannels, channel_names, build_pyramid!, images, labels,
    # Tables
    SpatialTable,
    tables, nobs, nvar, var_names, feature,
    # I/O
    SpatialDataZarr

include("coordsystems.jl")
include("dataset.jl")
include("elements.jl")
include("images.jl")
include("tables.jl")
include("views.jl")
include("zarr_io.jl")
include("show.jl")


end
