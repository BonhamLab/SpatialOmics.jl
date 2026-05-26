# ── Internal helper ────────────────────────────────────────────────────────────

_element_name(el::Union{SpatialPoints, SpatialShapes}) =
    (att = _dataset_ref(el); att === nothing ? "" : att[2])
_element_name(::Any) = ""

# ── Default dispatch — token inferred from argument types ─────────────────────

"""
    analyze(kind, args...) → SpatialRelation
    analyze(pts::SpatialPoints, cells::SpatialShapes)      → Expression relation
    analyze(src::SpatialShapes, dst::SpatialShapes)        → Membership relation
    analyze(rel::SpatialRelation{Expression}; k=30)        → KNN relation

Compute a spatial relation between elements.

Dispatch on the `RelationKind` token selects the algorithm:

- `analyze(Expression(), pts, cells)` — count transcripts per gene per cell.
  Each transcript is assigned to the first containing cell (bounding-box
  pre-filter, then exact point-in-polygon). Returns an n_cells × n_genes
  count matrix.
- `analyze(Membership(), src, dst)` — assign each point or shape in `src`
  to the containing shape in `dst`. `strict=true` requires full containment;
  default uses centroid or point containment.
- `analyze(KNN(k), rel)` — k-nearest-neighbour graph on an `Expression`
  relation. Requires `using NearestNeighbors`.

The two-argument forms (`pts, cells` and `src, dst`) infer the kind from
argument types and call the explicit form.

# See also
[`SpatialRelation`](@ref), [`Expression`](@ref), [`Membership`](@ref), [`KNN`](@ref), [`distances`](@ref)
"""
analyze(pts::SpatialPoints, cells::SpatialShapes) = analyze(Expression(), pts, cells)
analyze(src::SpatialShapes, dst::SpatialShapes)   = analyze(Membership(), src, dst)
analyze(rel::SpatialRelation{Expression}; k::Int=30) = analyze(KNN(k), rel)

# ── Expression — count transcripts per gene per cell ─────────────────────────
# Extent pre-filter (computed once) then GeometryOps.contains, same pattern as _mask.

function analyze(::Expression, pts::SpatialPoints, cells::SpatialShapes)
    n_cells = length(cells)
    n_genes = length(pts.feature_codebook)
    weights = zeros(Float32, n_cells, n_genes)
    extents = [_extent_of(cells.geometries[j]) for j in 1:n_cells]

    for i in eachindex(pts.coords)
        gid = pts.feature_id[i]
        gid == 0 && continue
        pt = pts.coords[i]
        px, py = Float64(pt[1]), Float64(pt[2])
        for j in 1:n_cells
            e = extents[j]
            e.xmin <= px <= e.xmax && e.ymin <= py <= e.ymax || continue
            GeometryOps.contains(cells.geometries[j], pt) || continue
            weights[j, gid] += 1f0
            break
        end
    end

    var_nt = isempty(pts.feature_codebook) ? NamedTuple() :
             NamedTuple{(:name,)}((pts.feature_codebook,))
    SpatialRelation(Expression(), _element_name(cells),
                    cells.instance_id, weights; var=var_nt)
end

# ── Membership — assign each source point/shape to a containing destination ───

function analyze(::Membership{strict}, pts::SpatialPoints, dst::SpatialShapes) where strict
    n_pts   = length(pts)
    n_dst   = length(dst)
    src_ids = Vector{Int32}(undef, n_pts)
    dst_ids = zeros(Int32, n_pts)
    extents = [_extent_of(dst.geometries[j]) for j in 1:n_dst]

    for i in 1:n_pts
        src_ids[i] = Int32(i)
        pt = pts.coords[i]
        px, py = Float64(pt[1]), Float64(pt[2])
        for j in 1:n_dst
            e = extents[j]
            e.xmin <= px <= e.xmax && e.ymin <= py <= e.ymax || continue
            GeometryOps.contains(dst.geometries[j], pt) || continue
            dst_ids[i] = dst.instance_id[j]
            break
        end
    end

    SpatialRelation(Membership{strict}(), _element_name(pts), _element_name(dst),
                    src_ids, dst_ids, nothing)
end

function analyze(::Membership{strict}, src::SpatialShapes, dst::SpatialShapes) where strict
    n_src   = length(src)
    n_dst   = length(dst)
    src_ids = copy(src.instance_id)
    dst_ids = zeros(Int32, n_src)
    dst_extents = [_extent_of(dst.geometries[j]) for j in 1:n_dst]

    for i in 1:n_src
        ctr = GeometryOps.centroid(src.geometries[i])
        cx, cy = Float64(ctr[1]), Float64(ctr[2])
        for j in 1:n_dst
            e = dst_extents[j]
            e.xmin <= cx <= e.xmax && e.ymin <= cy <= e.ymax || continue
            test_geom = strict ? src.geometries[i] : ctr
            GeometryOps.contains(dst.geometries[j], test_geom) || continue
            dst_ids[i] = dst.instance_id[j]
            break
        end
    end

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
    n = length(shapes_a)
    m = length(shapes_b)
    out = Vector{Float32}(undef, n)
    for i in 1:n
        ctr = GeometryOps.centroid(shapes_a.geometries[i])
        d   = Inf32
        for j in 1:m
            d = min(d, Float32(abs(GeometryOps.signed_distance(shapes_b.geometries[j], ctr))))
        end
        out[i] = d
    end
    out
end

function distances(shapes::SpatialShapes, roi::SpatialROI) :: Vector{Float32}
    [Float32(abs(GeometryOps.signed_distance(roi.geometry,
                 GeometryOps.centroid(shapes.geometries[i]))))
     for i in eachindex(shapes.geometries)]
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
