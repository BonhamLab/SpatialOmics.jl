module SpatialOmics

using GeometryBasics: Point, Point2f, Polygon, MultiPolygon, AbstractGeometry
using GeoInterface
using Tables
using StaticArrays
using OrderedCollections
using GeometryOps
using FlexiJoins
using Zarr
using JSON
using Parquet2
using CSV
using ImageBase: restrict
import ImageBase: scaleminmax
using MappedArrays: mappedarray
using Colors: Colorant, Gray, RGB
using Random: randperm
using SparseArrays
using CodecZstd
import ImageBase.ImageCore: colorview
using Images: load

export
    # Coordinate systems
    CoordinateSystem,
    # Transformations
    AbstractTransformation,
    Identity, Affine, Sequence,
    # translation, scaling, rotation, flip_y, compose — not exported; clash with Makie/LinearAlgebra.
    # Use SpatialOmics.translation(...) etc. when constructing pixel_to_cs transforms.
    apply, apply!, resolve,
    # Dataset
    BackingStore, SpatialDataset, AcquisitionSource,
    elements, coord_systems, transform, sources, source, source_attributes,
    with_dataset, keep!, save!, discard!, edit!, touch!, isdirty, dirty,
    # Elements
    SpatialPoints, SpatialShapes, SpatialShape,
    Polygon, MultiPolygon, Point2f,
    points, shapes,
    coords, features, feature_ids, origins, origin_ids, coord_system,
    geometries, instance_id, instance_ids, with_instance_ids,
    subsample, top_features, count_per_instance,
    # Views
    SpatialExtent, SpatialROI,
    SpatialElementView, SpatialDatasetView,
    geometry, roi, roi!,
    # Images + Labels
    SpatialImage, SpatialLabels, SpatialRasterTiles,
    data, nchannels, channel_names, build_pyramid!, ensure_pyramid!, images, labels,
    SpatialImageColorView, channel, scaleminmax, colorview, Gray, RGB,
    # Relations
    RelationKind, Membership, Expression,
    SpatialRelation,
    relations, source_ids, destination_ids, nobs, nvar, obs_names, var_names,
    annotate,
    # Analysis
    analyze, distances,
    PointDensity, density, ShapeColorView,
    # I/O
    SpatialDataZarr, native_store_version, write!,
    CosMx

include("coordsystems.jl")
include("relations.jl")
include("dataset.jl")
include("elements.jl")
include("views.jl")
include("images.jl")
include("analysis.jl")
include("zarr_io.jl")
include("show.jl")


end
