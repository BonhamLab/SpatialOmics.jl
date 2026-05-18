# ── RelationKind hierarchy ────────────────────────────────────────────────────
# Each kind is a singleton type so dispatch works: thing(rel) = thing(rel, rel.kind)

abstract type RelationKind end

# Membership{strict}: point/shape assigned to a containing shape
#   strict=false (default) — centroid-in-polygon or any overlap
#   strict=true            — full geometric containment
struct Membership{strict} <: RelationKind end
Membership(; strict::Bool=false) = Membership{strict}()

# Proximity: continuous pairwise distance matrix (symmetric, n×n)
# Boolean within-r case: apply weights .< r after the fact
struct Proximity <: RelationKind end

# KNN: k-nearest-neighbor graph between elements of the same collection
# weights :: Matrix{Float32} n×k, dst_ids is flat length n*k edge list
struct KNN <: RelationKind
    k :: Int
end
KNN(; k::Int = 30) = KNN(k)

# Expression: bipartite weighted relation (cells × genes count matrix)
# replaces the former SpatialTable type
struct Expression <: RelationKind end

# ── SpatialRelation ───────────────────────────────────────────────────────────

struct SpatialRelation{K<:RelationKind, W}
    src     :: String                    # source element name in ds.elements
    dst     :: Union{String, Nothing}    # dest element name; Nothing for Expression
    src_ids :: Vector{Int32}             # instance_ids of source rows/nodes
    dst_ids :: Vector{Int32}             # instance_ids of dest (empty for Expression)
    weights :: W                         # Nothing | AbstractMatrix
    obs     :: Any                       # per-source-row metadata (Tables-compatible)
    var     :: Any                       # per-dest-col metadata (Expression only)
    kind    :: K                         # singleton enables dispatch
end

# Inner constructor — infers K and W from kind and weights
function SpatialRelation(src::String, dst, src_ids, dst_ids, weights, obs, var, kind::K) where {K<:RelationKind}
    SpatialRelation{K, typeof(weights)}(
        src, dst,
        Vector{Int32}(src_ids),
        Vector{Int32}(dst_ids),
        weights, obs, var, kind)
end

# ── Convenience constructors ──────────────────────────────────────────────────

# Membership / Proximity / KNN: no var metadata
function SpatialRelation(kind::RelationKind, src::String, dst::String,
                          src_ids, dst_ids, weights=nothing;
                          obs=NamedTuple())
    SpatialRelation(src, dst, src_ids, dst_ids, weights, obs, NamedTuple(), kind)
end

# Expression: has obs (per-cell) and var (per-gene)
function SpatialRelation(::Expression, src::String,
                          src_ids, weights::AbstractMatrix;
                          obs=NamedTuple(), var=NamedTuple())
    SpatialRelation(src, nothing, src_ids, Int32[], weights, obs, var, Expression())
end

# ── Accessors ─────────────────────────────────────────────────────────────────

nobs(rel::SpatialRelation) = length(rel.src_ids)

function nvar(rel::SpatialRelation{Expression})
    isnothing(rel.weights) ? 0 : size(rel.weights, 2)
end

function var_names(rel::SpatialRelation{Expression})
    t = Tables.columns(rel.var)
    Tables.columnnames(t) !== () && hasproperty(t, :name) ?
        collect(t.name) : string.(1:nvar(rel))
end

# ── annotate — pure: returns new SpatialRelation with obs column added ────────

function annotate(rel::SpatialRelation, labels::AbstractVector; key::Symbol)
    length(labels) == nobs(rel) ||
        error("labels length ($(length(labels))) ≠ nobs(rel) ($(nobs(rel)))")
    new_obs = merge(Tables.columntable(rel.obs), NamedTuple{(key,)}((labels,)))
    SpatialRelation(rel.src, rel.dst, rel.src_ids, rel.dst_ids,
                    rel.weights, new_obs, rel.var, rel.kind)
end
