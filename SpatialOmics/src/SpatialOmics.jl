module SpatialOmics

using GeometryBasics: Point, Point2f, Polygon, AbstractGeometry
using GeoInterface
using Tables
using StaticArrays
using OrderedCollections
using GeometryOps
using Zarr
using JSON
using Parquet2
using CSV
using ImageBase: restrict
import ImageBase: scaleminmax
using Colors: Colorant
using Random: randperm
import ImageBase.ImageCore: colorview

export
    # Coordinate systems
    CoordinateSystem,
    # Transformations
    AbstractTransformation,
    Identity, Affine, Sequence,
    translation, scaling, rotation, flip_y, compose,
    apply, apply!, resolve,
    # Dataset
    BackingStore, SpatialDataset,
    elements, coord_systems, transform,
    with_dataset, keep!,
    # Elements
    SpatialPoints, SpatialShapes,
    points, shapes,
    coords, features, feature_ids, coord_system,
    geometries, bbox, instance_id, instance_ids,
    subsample, top_features, count_per_instance,
    # Views
    SpatialExtent, SpatialROI,
    SpatialElementView, SpatialDatasetView,
    geometry,
    # Images + Labels
    SpatialImage, SpatialLabels,
    data, nchannels, channel_names, build_pyramid!, images, labels,
    SpatialImageColorView, channel, scaleminmax, colorview,
    # Tables
    SpatialTable,
    tables, nobs, nvar, var_names, feature,
    # I/O
    SpatialDataZarr,
    CosMx

include("coordsystems.jl")
include("dataset.jl")
include("elements.jl")
include("views.jl")
include("images.jl")
include("tables.jl")
include("zarr_io.jl")
include("show.jl")


end
