# ---------------------------------------------------------------------------
# statistics/SpatialStatistics.jl — included directly into SpatialAnalysis
# ---------------------------------------------------------------------------

using NearestNeighbors
using GeometryBasics
using GeometryOps
using LinearAlgebra

import SpatialOmicsBase: _geom_centroid

# ---------------------------------------------------------------------------
# distances — distance from points/cells to a region or border
# ---------------------------------------------------------------------------

"""
    distances(from::SpatialPoints, to::SpatialShapes; signed::Bool=false) -> Vector{Float64}
    distances(from::SpatialShapes, to::SpatialShapes; signed::Bool=false) -> Vector{Float64}
    distances(from::SpatialPoints, to::SpatialExtent)                     -> Vector{Float64}
    distances(from::SpatialShapes, to::SpatialExtent)                     -> Vector{Float64}

Compute the distance from each point or cell centroid in `from` to the
nearest boundary in `to`.

When `to` is a `SpatialShapes` element, a KDTree is built over the shape
centroids (O(N log M)) and the exact boundary distance is computed via
`GeometryOps` against the nearest shape only. Use `signed=true` to get
negative distances for points/cells inside a polygon.

When `to` is a `SpatialExtent`, returns the distance to the nearest edge of
the bounding box (interior points return their distance to the nearest edge,
not zero).

`from::SpatialShapes` uses the `x_centroid`/`y_centroid` columns in
`from.features` as source coordinates.

Returns a `Vector{Float64}` aligned to the rows of `from`.
"""
function distances(from::SpatialPoints, to::SpatialShapes; signed::Bool=false)
    n_pts = size(from.coordinates, 1)
    n_shapes = length(to.geometries)
    n_shapes == 0 && return fill(Inf, n_pts)

    tree, cx, cy = _centroid_tree(to)
    dist_fn = signed ? GeometryOps.signed_distance : GeometryOps.distance

    result = ThreadsX.map(1:n_pts) do i
        x = Float64(from.coordinates[i, 1])
        y = Float64(from.coordinates[i, 2])
        idx = _nearest_shape(tree, x, y)
        _shape_distance(to.geometries[idx], dist_fn, x, y, cx[idx], cy[idx])
    end
    return result
end

function distances(from::SpatialShapes, to::SpatialShapes; signed::Bool=false)
    feats = from.features
    n = length(from.geometries)
    n == 0 && return Float64[]
    n_to = length(to.geometries)
    n_to == 0 && return fill(Inf, n)

    has_c = ncol(feats) > 0 &&
            "x_centroid" in names(feats) &&
            "y_centroid" in names(feats)

    tree, cx, cy = _centroid_tree(to)
    dist_fn = signed ? GeometryOps.signed_distance : GeometryOps.distance

    result = ThreadsX.map(1:n) do i
        if has_c
            x, y = Float64(feats.x_centroid[i]), Float64(feats.y_centroid[i])
        else
            c = _geom_centroid(from.geometries[i])
            c === nothing && return Inf
            x, y = Float64(c[1]), Float64(c[2])
        end
        idx = _nearest_shape(tree, x, y)
        _shape_distance(to.geometries[idx], dist_fn, x, y, cx[idx], cy[idx])
    end
    return result
end

function distances(from::SpatialPoints, to::SpatialExtent)
    n = size(from.coordinates, 1)
    return [_dist_to_extent(Float64(from.coordinates[i,1]),
                            Float64(from.coordinates[i,2]), to) for i in 1:n]
end

function distances(from::SpatialShapes, to::SpatialExtent)
    feats = from.features
    n = length(from.geometries)
    n == 0 && return Float64[]
    has_c = ncol(feats) > 0 &&
            "x_centroid" in names(feats) &&
            "y_centroid" in names(feats)
    return map(1:n) do i
        if has_c
            x, y = Float64(feats.x_centroid[i]), Float64(feats.y_centroid[i])
        else
            c = _geom_centroid(from.geometries[i])
            c === nothing && return Inf
            x, y = Float64(c[1]), Float64(c[2])
        end
        _dist_to_extent(x, y, to)
    end
end

# ---------------------------------------------------------------------------
# distances internals
# ---------------------------------------------------------------------------

# Build a KDTree over the centroids of `shp`. Returns (tree, cx, cy).
function _centroid_tree(shp::SpatialShapes)
    n = length(shp.geometries)
    feats = shp.features
    has_c = ncol(feats) > 0 &&
            "x_centroid" in names(feats) &&
            "y_centroid" in names(feats)
    cx = Vector{Float64}(undef, n)
    cy = Vector{Float64}(undef, n)
    if has_c
        cx .= Float64.(feats.x_centroid)
        cy .= Float64.(feats.y_centroid)
    else
        for i in 1:n
            c = _geom_centroid(shp.geometries[i])
            cx[i] = c !== nothing ? Float64(c[1]) : 0.0
            cy[i] = c !== nothing ? Float64(c[2]) : 0.0
        end
    end
    data = Matrix{Float64}(undef, 2, n)
    data[1, :] .= cx
    data[2, :] .= cy
    return KDTree(data), cx, cy
end

# Find index of nearest shape given a KDTree built from _centroid_tree.
function _nearest_shape(tree::KDTree, x::Float64, y::Float64)
    idxs, _ = knn(tree, [x, y], 1)
    return idxs[1]
end

# Compute exact distance from point (x,y) to a geometry using dist_fn.
# Falls back to centroid distance for unsupported geometry types.
function _shape_distance(geom, dist_fn, x::Float64, y::Float64, cx::Float64, cy::Float64)
    pt = GeometryBasics.Point2(x, y)
    if geom isa GeometryBasics.Circle
        # Analytic: signed distance to circle boundary
        d = hypot(cx - x, cy - y) - geom.r
        return dist_fn === GeometryOps.signed_distance ? d : max(d, 0.0)
    end
    return dist_fn(pt, geom)
end

# Distance from point (x,y) to nearest edge of the bounding box.
# Points inside return their distance to the nearest wall; points outside
# return the Euclidean distance to the nearest corner/edge.
function _dist_to_extent(x::Float64, y::Float64, ext::SpatialExtent)
    dx = max(Float64(ext.xmin) - x, 0.0, x - Float64(ext.xmax))
    dy = max(Float64(ext.ymin) - y, 0.0, y - Float64(ext.ymax))
    if dx == 0.0 && dy == 0.0
        # inside: distance to nearest edge
        return min(x - Float64(ext.xmin), Float64(ext.xmax) - x,
                   y - Float64(ext.ymin), Float64(ext.ymax) - y)
    end
    return hypot(dx, dy)
end

# ---------------------------------------------------------------------------
# classify_cells — assign cell-type labels from markers or signatures
# ---------------------------------------------------------------------------

"""
    classify_cells(tbl::SpatialTable, markers::Dict{String,Vector{String}};
                   method::Symbol = :mean_expression,
                   min_markers::Int = 1,
                   var_gene_key::String = first(names(tbl.var))) -> Vector{String}

    classify_cells(tbl::SpatialTable, signatures::AbstractMatrix;
                   var_names::Vector{String},
                   type_names::Vector{String},
                   method::Symbol = :dot_product,
                   var_gene_key::String = first(names(tbl.var))) -> Vector{String}

Assign a cell-type label to each observation in `tbl`.

**Marker-based** (`markers` form): each cell is scored by the mean expression
of the marker genes for each type; the type with the highest score is assigned.
Cells with a max score of zero are labelled `"unknown"`. Supported `method`:
`:mean_expression`.

**Signature-based** (`signatures` form): scores are computed as the dot product
(or Pearson correlation) of each cell's expression vector against each type's
reference signature. Supported `method` values: `:dot_product`, `:correlation`.

`min_markers` (marker form only): minimum number of marker genes that must be
present in the table for a cell type to be scored; types with fewer matches
receive a score of 0.

`var_gene_key`: column of `tbl.var` containing gene/feature names
(default: first column).

Returns a `Vector{String}` of length `nrow(tbl.obs)`. Does not mutate `tbl`.
Add the result to `tbl.obs` as needed:
```julia
tbl.obs[!, "cell_type"] = classify_cells(tbl, markers)
```
"""
function classify_cells(
    tbl::SpatialTable,
    markers::Dict{String,Vector{String}};
    method::Symbol = :mean_expression,
    min_markers::Int = 1,
    var_gene_key::String = first(names(tbl.var)),
)
    gene_names = tbl.var[!, var_gene_key]
    gene_index = Dict{String,Int}(g => i for (i, g) in enumerate(gene_names))

    type_names = collect(keys(markers))
    n_types = length(type_names)
    n_obs   = size(tbl.data, 1)
    X = Matrix{Float64}(tbl.data)
    scores = zeros(Float64, n_obs, n_types)

    for (j, t) in enumerate(type_names)
        gene_idxs = [gene_index[g] for g in markers[t] if haskey(gene_index, g)]
        length(gene_idxs) < min_markers && continue
        col_sum = sum(view(X, :, gene_idxs), dims=2)
        scores[:, j] .= vec(col_sum) ./ length(gene_idxs)
    end

    labels = Vector{String}(undef, n_obs)
    for i in 1:n_obs
        best = argmax(view(scores, i, :))
        labels[i] = scores[i, best] > 0 ? type_names[best] : "unknown"
    end
    return labels
end

function classify_cells(
    tbl::SpatialTable,
    signatures::AbstractMatrix;
    var_names::Vector{String},
    type_names::Vector{String},
    method::Symbol = :dot_product,
    var_gene_key::String = first(names(tbl.var)),
)
    size(signatures, 1) == length(var_names) ||
        throw(ArgumentError("signatures rows ($(size(signatures,1))) must match length(var_names) ($(length(var_names)))"))
    size(signatures, 2) == length(type_names) ||
        throw(ArgumentError("signatures columns ($(size(signatures,2))) must match length(type_names) ($(length(type_names)))"))

    tbl_genes = tbl.var[!, var_gene_key]
    sig_index = Dict{String,Int}(g => i for (i, g) in enumerate(var_names))

    shared     = [g for g in tbl_genes if haskey(sig_index, g)]
    tbl_cols   = [findfirst(==(g), tbl_genes) for g in shared]
    sig_rows   = [sig_index[g] for g in shared]

    X_sub = Matrix{Float64}(tbl.data[:, tbl_cols])
    S_sub = Matrix{Float64}(signatures[sig_rows, :])

    if method == :dot_product
        score_mat = X_sub * S_sub
    elseif method == :correlation
        X_z = (X_sub .- mean(X_sub, dims=1)) ./ (std(X_sub, dims=1) .+ 1e-8)
        S_z = (S_sub .- mean(S_sub, dims=1)) ./ (std(S_sub, dims=1) .+ 1e-8)
        score_mat = X_z * S_z
    else
        throw(ArgumentError("classify_cells: unknown method :$method (choose :dot_product or :correlation)"))
    end

    n_obs = size(score_mat, 1)
    return [type_names[argmax(score_mat[i, :])] for i in 1:n_obs]
end

"""
    spatial_autocorrelation(ds::SpatialDataset,
                            feature::String;
                            method::Symbol = :morans_i,
                            n_permutations::Int = 999) -> NamedTuple

Compute a spatial autocorrelation statistic for `feature` across the spatial
observations in `ds`.

Returns a NamedTuple with fields `statistic`, `pvalue`, `zscore`.

Supported `method` values (planned):
- `:morans_i`   — Moran's I global autocorrelation
- `:gearys_c`   — Geary's C global autocorrelation
- `:getis_ord`  — Getis-Ord G* local statistic (hotspot variant)
"""
function spatial_autocorrelation(
    ds::SpatialDataset,
    feature::String;
    method::Symbol = :morans_i,
    n_permutations::Int = 999,
)
    # TODO Phase 2: implement; parallelise permutations with ThreadsX
    error("spatial_autocorrelation: not yet implemented (Phase 2 deliverable)")
end

"""
    hotspot_detection(ds::SpatialDataset,
                      feature::String;
                      fdr_threshold::Float64 = 0.05) -> DataFrame

Identify statistically significant spatial hotspots (high-expression regions)
for `feature`. Returns a per-observation DataFrame with columns:
`local_statistic`, `pvalue`, `fdr`, `is_hotspot`.
"""
function hotspot_detection(
    ds::SpatialDataset,
    feature::String;
    fdr_threshold::Float64 = 0.05,
)
    # TODO Phase 2: implement
    error("hotspot_detection: not yet implemented (Phase 2 deliverable)")
end

"""
    spatial_clustering(ds::SpatialDataset;
                       method::Symbol = :bayesspace,
                       n_clusters::Union{Int, Nothing} = nothing,
                       use_pca::Bool = true,
                       n_pcs::Int = 30) -> Vector{Int}

Cluster spatial observations into coherent tissue domains.
Returns a vector of integer cluster labels, one per observation.

Supported `method` values (planned):
- `:bayesspace`   — BayesSpace spatially-aware clustering (via PythonCall)
- `:leiden`       — Leiden graph clustering on spatial neighbourhood graph
- `:kmeans`       — k-means on PCA-reduced expression + spatial coordinates
"""
function spatial_clustering(
    ds::SpatialDataset;
    method::Symbol = :bayesspace,
    n_clusters::Union{Int, Nothing} = nothing,
    use_pca::Bool = true,
    n_pcs::Int = 30,
)
    # TODO Phase 2: implement
    error("spatial_clustering: not yet implemented (Phase 2 deliverable)")
end

