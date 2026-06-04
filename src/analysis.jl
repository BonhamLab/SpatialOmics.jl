# ── Internal helper ────────────────────────────────────────────────────────────

_element_name(el::Union{SpatialPoints, SpatialShapes}) =
    (att = _dataset_ref(el); att === nothing ? "" : att[2])
_element_name(::Any) = ""

# ── Default dispatch — token inferred from argument types ─────────────────────

"""
    analyze(kind, args...) → SpatialRelation
    analyze(pts::SpatialPoints, cells::SpatialShapes)  → Expression relation
    analyze(src::SpatialShapes, dst::SpatialShapes)    → Membership relation

Compute a spatial relation between elements.

Dispatch on the `RelationKind` token selects the algorithm:

- `analyze(Expression(), pts, cells)` — count transcripts per gene per cell.
  Each transcript is assigned to the first containing cell (bounding-box
  pre-filter, then exact point-in-polygon). Returns an n_cells × n_genes
  count matrix.
- `analyze(Membership(), src, dst)` — assign each point or shape in `src`
  to the containing shape in `dst`. `strict=true` requires full containment;
  default uses centroid or point containment.

The two-argument forms (`pts, cells` and `src, dst`) infer the kind from
argument types and call the explicit form.

# See also
[`SpatialRelation`](@ref), [`Expression`](@ref), [`Membership`](@ref), [`distances`](@ref)
"""
analyze(pts::SpatialPoints, cells::SpatialShapes; kw...) = analyze(Expression(), pts, cells; kw...)
analyze(src::SpatialShapes, dst::SpatialShapes; kw...)   = analyze(Membership(), src, dst; kw...)

# ── Multi-level: aggregate an obs column across a spatial grouping ────────────
# analyze(src_shapes, dst_shapes; obs, by) → Expression relation where rows are
# dst groups and columns are the unique values of obs[by].
# Example: cells → ROIs grouped by :cell_type gives an n_ROIs × n_types matrix.

"""
    analyze(src, dst; obs, by, predicate=GO.within) → SpatialRelation{Expression}

Second-level aggregation: group `src` shapes into `dst` shapes, then count by
the `by` column of `obs`.

```julia
rel_rois = analyze(cells, rois; obs=rel_cells.obs, by=:cell_type)
rel_rois[:, "Epithelial"]   # count of Epithelial cells per ROI
```

# See also
[`analyze`](@ref), [`annotate`](@ref)
"""
function analyze(src::SpatialShapes, dst::SpatialShapes, obs;
                 by::Symbol, predicate=GeometryOps.within)
    labels    = getproperty(Tables.columntable(obs), by)
    uniq      = sort(unique(labels))
    label_pos = Dict(l => i for (i, l) in enumerate(uniq))
    n_dst     = length(dst)
    n_labels  = length(uniq)
    weights   = zeros(Float32, n_dst, n_labels)
    dst_pos   = Dict{Int32,Int}(id => i for (i, id) in enumerate(dst.instance_id))

    src_geoms = GeometryOps.centroid.(src.geometries)
    src_tbl   = Tables.rowtable((geom=src_geoms,         label=labels,
                                  src_instance_id=src.instance_id))
    dst_tbl   = Tables.rowtable((poly=dst.geometries,    dst_instance_id=dst.instance_id))

    for row in innerjoin((src_tbl, dst_tbl), by_pred(:geom, predicate, :poly))
        lpos = get(label_pos, row[1].label, 0)
        lpos == 0 && continue
        dpos = get(dst_pos, row[2].dst_instance_id, 0)
        dpos == 0 && continue
        weights[dpos, lpos] += 1f0
    end

    var_nt = NamedTuple{(:name,)}((collect(String, string.(uniq)),))
    SpatialRelation(Expression(), _element_name(dst),
                    dst.instance_id, weights; var=var_nt)
end

# ── Expression — count transcripts per gene per cell ─────────────────────────
# FlexiJoins does the spatial join (with tree-indexed polygons); we accumulate
# into the count matrix in a single pass over the matched pairs.
# Points go first (simpler geom), shapes second (tree-indexed).
# GO.within(point, polygon) = "point is within polygon".

function analyze(::Expression, pts::SpatialPoints, cells::SpatialShapes;
                 predicate=GeometryOps.within)
    n_cells  = length(cells)
    n_genes  = length(pts.feature_codebook)
    weights  = zeros(Float32, n_cells, n_genes)
    cell_pos = Dict{Int32,Int}(id => i for (i, id) in enumerate(cells.instance_id))

    # FlexiJoins expects row-iterable tables; result rows are Tuple{src_row, dst_row}.
    # Points go first (simpler geoms), cells second (tree-indexed by FlexiJoins).
    pts_tbl   = Tables.rowtable((pt=pts.coords,         feature_id=pts.feature_id))
    cells_tbl = Tables.rowtable((poly=cells.geometries, instance_id=cells.instance_id))

    for row in innerjoin((pts_tbl, cells_tbl), by_pred(:pt, predicate, :poly))
        gid  = row[1].feature_id
        gid == 0 && continue
        cpos = get(cell_pos, row[2].instance_id, 0)
        cpos == 0 && continue
        weights[cpos, gid] += 1f0
    end

    var_nt = isempty(pts.feature_codebook) ? NamedTuple() :
             NamedTuple{(:name,)}((pts.feature_codebook,))
    SpatialRelation(Expression(), _element_name(cells),
                    cells.instance_id, weights; var=var_nt)
end

# ── Membership — assign each source point/shape to a containing destination ───
# strict=false (default): point-in-polygon for transcripts, centroid-in-polygon
# for shapes. strict=true: full geometric containment of src within dst shape.
# Returns only matched pairs (innerjoin semantics — no 0-sentinel rows).

function analyze(::Membership{strict}, pts::SpatialPoints, dst::SpatialShapes;
                 predicate=GeometryOps.within) where strict
    pts_tbl = Tables.rowtable((pt=pts.coords,          pos=Int32.(eachindex(pts.coords))))
    dst_tbl = Tables.rowtable((poly=dst.geometries,    dst_instance_id=dst.instance_id))

    joined  = collect(innerjoin((pts_tbl, dst_tbl), by_pred(:pt, predicate, :poly)))
    src_ids = Int32[row[1].pos            for row in joined]
    dst_ids = Int32[row[2].dst_instance_id for row in joined]

    SpatialRelation(Membership{strict}(), _element_name(pts), _element_name(dst),
                    src_ids, dst_ids, nothing)
end

function analyze(::Membership{strict}, src::SpatialShapes, dst::SpatialShapes;
                 predicate=nothing) where strict
    # strict=false: centroid-in-polygon; strict=true: full shape containment.
    src_geoms = strict ? src.geometries : GeometryOps.centroid.(src.geometries)
    pred      = isnothing(predicate) ? GeometryOps.within : predicate

    src_tbl = Tables.rowtable((geom=src_geoms,          src_instance_id=src.instance_id))
    dst_tbl = Tables.rowtable((poly=dst.geometries,      dst_instance_id=dst.instance_id))

    joined  = collect(innerjoin((src_tbl, dst_tbl), by_pred(:geom, pred, :poly)))
    src_ids = Int32[row[1].src_instance_id for row in joined]
    dst_ids = Int32[row[2].dst_instance_id for row in joined]

    SpatialRelation(Membership{strict}(), _element_name(src), _element_name(dst),
                    src_ids, dst_ids, nothing)
end

# ── distances ─────────────────────────────────────────────────────────────────

"""
    distances(shapes_a, shapes_b) → Vector{Float32}
    distances(shapes, roi) → Vector{Float32}

Compute the minimum distance from each shape in `shapes_a` to the nearest shape
in `shapes_b` (or to a `SpatialROI` boundary), using centroid-to-shape
signed distance.

Returns a `Float32` vector of length `length(shapes_a)`.

# See also
[`analyze`](@ref), [`Proximity`](@ref)
"""
function distances(shapes_a::SpatialShapes, shapes_b::SpatialShapes) :: Vector{Float32}
    [Float32(minimum(GeometryOps.distance(GeometryOps.centroid(g_a), g_b)
                     for g_b in geometries(shapes_b)))
     for g_a in geometries(shapes_a)]
end

function distances(shapes::SpatialShapes, roi::SpatialROI) :: Vector{Float32}
    # signed_distance gives negative values inside — abs gives boundary distance
    # even for shapes whose centroid is inside the ROI.
    [Float32(abs(GeometryOps.signed_distance(roi.geometry,
                 GeometryOps.centroid(g))))
     for g in geometries(shapes)]
end

# ── PointDensity ──────────────────────────────────────────────────────────────

"""
    PointDensity

Lazy density estimate descriptor for a `SpatialPoints` collection.

Produced by `density(pts; resolution, feature)`. Passed to Makie plot verbs to
render a rasterised kernel density map at display time.

# See also
[`density`](@ref)
"""
struct PointDensity
    pts        :: SpatialPoints
    resolution :: Int
    feature    :: Union{String, Nothing}
end

"""
    density(pts; resolution=512, feature=nothing) → PointDensity

Create a lazy density estimate descriptor for `pts`.

`resolution` sets the output grid size (pixels along the longer axis).
`feature` restricts density computation to a single feature label; `nothing`
uses all points.

# See also
[`PointDensity`](@ref)
"""
density(pts::SpatialPoints; resolution::Int=512, feature::Union{String,Nothing}=nothing) =
    PointDensity(pts, resolution, feature)

# ── ShapeColorView ────────────────────────────────────────────────────────────

"""
    ShapeColorView

Display descriptor that pairs a `SpatialShapes` collection with per-shape color
values from a `SpatialRelation`.

Passed to Makie's `poly!` to render shapes coloured by an expression or other
quantitative measure.

# See also
[`SpatialShapes`](@ref), [`SpatialRelation`](@ref)
"""
struct ShapeColorView
    shapes   :: SpatialShapes
    rel      :: SpatialRelation
    color_by :: Symbol
    colormap :: Any
end
