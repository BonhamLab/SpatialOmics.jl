"""
    Preprocessing

Count normalisation, batch effect correction, missing-value imputation,
and spatial cropping for spatial transcriptomics datasets. Functions are
designed to be pure and composable, operating on `SpatialDataset` objects.
"""
module Preprocessing

using DataFrames
using Statistics
using GeometryBasics
using GeometryOps
using Tables
using SpatialOmicsBase

import SpatialOmicsBase: _crop, _geom_centroid, crop

export normalize_counts, correct_batch_effects, impute_missing

# ---------------------------------------------------------------------------
# crop — materialise a spatially-filtered copy of an element or dataset
# ---------------------------------------------------------------------------

"""
    crop(pts::SpatialPoints, ext::SpatialExtent) -> SpatialPoints
    crop(shp::SpatialShapes, ext::SpatialExtent) -> SpatialShapes
    crop(tbl::SpatialTable,  ext::SpatialExtent, shp::SpatialShapes) -> SpatialTable
    crop(ds::SpatialDataset, ext::SpatialExtent; ...) -> SpatialDataset

    crop(pts::SpatialPoints, roi::SpatialShapes) -> SpatialPoints
    crop(shp::SpatialShapes, roi::SpatialShapes) -> SpatialShapes

Materialise a spatially-filtered copy of the element(s) within `ext` (an
axis-aligned bounding box) or `roi` (arbitrary polygon region — union of all
geometries in the shapes element).

This is the explicit materialisation path. For lazy, non-copying views use
`view(element, ext)` instead.

# Dataset keyword arguments
- `points_keys`  : keys to filter; `nothing` (default) filters all
- `shapes_keys`  : keys to filter; `nothing` (default) filters all
- `tables_keys`  : keys to filter; `nothing` (default) filters all

Images and labels are left untouched (use lazy views for pixel cropping).
"""
crop(pts::SpatialPoints,  ext::SpatialExtent) = _crop(pts, ext)
crop(shp::SpatialShapes,  ext::SpatialExtent) = _crop(shp, ext)

function crop(tbl::SpatialTable, ext::SpatialExtent, shp::SpatialShapes)
    cropped_shp = _crop(shp, ext)
    ikey        = Symbol(instance_key(tbl))
    valid_ids   = Set(Tables.getcolumn(cropped_shp, ikey))
    mask = [id in valid_ids for id in tbl.obs[!, String(ikey)]]
    return SpatialTable(tbl.data[mask, :], tbl.obs[mask, :], tbl.var, tbl.metadata)
end

function crop(ds::SpatialDataset{T}, ext::SpatialExtent;
                points_keys=nothing, shapes_keys=nothing, tables_keys=nothing) where T
    pkeys = points_keys === nothing ? collect(keys(ds.points)) : points_keys
    skeys = shapes_keys  === nothing ? collect(keys(ds.shapes))  : shapes_keys
    tkeys = tables_keys  === nothing ? collect(keys(ds.tables))  : tables_keys

    new_points = Dict(k => _crop(ds.points[k], ext)  for k in pkeys)
    new_shapes = Dict(k => _crop(ds.shapes[k], ext)  for k in skeys)
    new_tables = Dict{String, SpatialTable}()
    for k in tkeys
        tbl = ds.tables[k]
        r = region(tbl)
        if r !== nothing && haskey(ds.shapes, r)
            new_tables[k] = crop(tbl, ext, ds.shapes[r])
        else
            new_tables[k] = tbl  # no linked shapes; include as-is
        end
    end

    return SpatialDataset{T}(
        ds.images,
        ds.labels,
        new_points,
        new_shapes,
        new_tables,
        ds.coordinate_systems,
        ds.transformations,
        ds.metadata,
    )
end

# Polygon ROI variants — test containment using GeometryOps.signed_distance.
# Multiple geometries in `roi` are treated as a union (inside any = included).

# Bounding box of a SpatialShapes element (union of all geometry bounds).
function _shapes_extent(roi::SpatialShapes)
    xs = Float64[]; ys = Float64[]
    for s in roi.shapes
        _collect_coords!(xs, ys, s.geometry)
    end
    isempty(xs) && return SpatialExtent(-Inf, Inf, -Inf, Inf)
    return SpatialExtent(minimum(xs), maximum(xs), minimum(ys), maximum(ys))
end

_collect_coords!(xs, ys, g::GeometryBasics.Circle) =
    (push!(xs, g.center[1] - g.r, g.center[1] + g.r);
     push!(ys, g.center[2] - g.r, g.center[2] + g.r))

function _collect_coords!(xs, ys, g::GeometryBasics.Polygon)
    for pt in GeometryBasics.coordinates(g)
        push!(xs, Float64(pt[1])); push!(ys, Float64(pt[2]))
    end
end

_collect_coords!(xs, ys, ::Any) = nothing  # unsupported geometry type — skip

function crop(pts::SpatialPoints, roi::SpatialShapes)
    # Pre-filter to extent bounding box to reduce per-point polygon tests.
    roi_ext = _shapes_extent(roi)
    candidates = _crop(pts, roi_ext)
    n = size(candidates.coordinates, 1)
    mask = Vector{Bool}(undef, n)
    for i in 1:n
        x, y = Float64(candidates.coordinates[i, 1]), Float64(candidates.coordinates[i, 2])
        pt = GeometryBasics.Point2(x, y)
        mask[i] = any(GeometryOps.signed_distance(pt, s.geometry) <= 0 for s in roi.shapes)
    end
    feats = candidates.features isa DataFrame ? candidates.features : DataFrame(candidates.features)
    return SpatialPoints(candidates.coordinates[mask, :], feats[mask, :], candidates.metadata)
end

function crop(shp::SpatialShapes, roi::SpatialShapes)
    roi_ext    = _shapes_extent(roi)
    candidates = _crop(shp, roi_ext)
    n_cands    = length(candidates.shapes)
    mask       = Vector{Bool}(undef, n_cands)
    for i in 1:n_cands
        cxy = SpatialOmicsBase._geom_centroid(candidates.shapes[i].geometry)
        cxy === nothing && (mask[i] = false; continue)
        pt = GeometryBasics.Point2(Float64(cxy[1]), Float64(cxy[2]))
        mask[i] = any(GeometryOps.signed_distance(pt, s.geometry) <= 0 for s in roi.shapes)
    end
    return SpatialShapes(candidates.shapes[mask], candidates.metadata)
end

# ---------------------------------------------------------------------------
# Stub functions (Phase 2 planned)
# ---------------------------------------------------------------------------

"""
    normalize_counts(ds::SpatialDataset;
                     method::Symbol = :scran,
                     target_sum::Real = 1e4) -> SpatialDataset

Normalise raw count data in `ds`.

Supported `method` values (planned):
- `:scran`    — pooling-based size-factor normalisation (via PythonCall bridge)
- `:total`    — divide each cell by its total count, multiply by `target_sum`
- `:log1p`    — log1p-transform after total normalisation
"""
function normalize_counts(
    ds::SpatialDataset;
    method::Symbol = :scran,
    target_sum::Real = 1e4,
)
    # TODO Phase 2: implement
    error("normalize_counts: not yet implemented (Phase 2 deliverable)")
end

"""
    correct_batch_effects(ds::SpatialDataset;
                          batch_key::String = "sample",
                          method::Symbol = :harmony) -> SpatialDataset

Remove technical batch variation across slides or samples.

Supported `method` values (planned):
- `:harmony`  — Harmony integration (via PythonCall)
- `:combat`   — ComBat parametric batch correction
- `:scvi`     — scVI deep generative model (via PythonCall)
"""
function correct_batch_effects(
    ds::SpatialDataset;
    batch_key::String = "sample",
    method::Symbol = :harmony,
)
    # TODO Phase 2: implement
    error("correct_batch_effects: not yet implemented (Phase 2 deliverable)")
end

"""
    impute_missing(ds::SpatialDataset;
                   method::Symbol = :tangram,
                   reference::Union{Nothing, SpatialDataset} = nothing)
        -> SpatialDataset

Impute missing gene expression values using spatial context.

Supported `method` values (planned):
- `:tangram`  — Tangram transcript-level mapping (requires reference)
- `:gimvi`    — GIMVI spatial imputation
- `:knn`      — simple k-NN smoothing over spatial neighbours
"""
function impute_missing(
    ds::SpatialDataset;
    method::Symbol = :tangram,
    reference::Union{Nothing, SpatialDataset} = nothing,
)
    # TODO Phase 2: implement
    error("impute_missing: not yet implemented (Phase 2 deliverable)")
end

end # module Preprocessing
