# ── Internal helper ────────────────────────────────────────────────────────────

_element_name(el::Union{SpatialPoints, SpatialShapes}) =
    (att = _dataset_ref(el); att === nothing ? "" : att[2])
_element_name(::Any) = ""

# ── Default dispatch — token inferred from argument types ─────────────────────

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

# ── KNN stub — actual implementation lives in ext/NearestNeighborsExt.jl ──────

function analyze(::KNN, ::SpatialRelation{Expression})
    error("KNN requires a backend. Load one:\n  using NearestNeighbors")
end

# ── distances ─────────────────────────────────────────────────────────────────

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

struct PointDensity
    pts        :: SpatialPoints
    resolution :: Int
    feature    :: Union{String, Nothing}
end

density(pts::SpatialPoints; resolution::Int=512, feature::Union{String,Nothing}=nothing) =
    PointDensity(pts, resolution, feature)

# ── ShapeColorView ────────────────────────────────────────────────────────────

struct ShapeColorView
    shapes   :: SpatialShapes
    rel      :: SpatialRelation
    color_by :: Symbol
    colormap :: Any
end
