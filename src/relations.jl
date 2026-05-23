# ── RelationKind hierarchy ────────────────────────────────────────────────────

"""
    RelationKind

Abstract supertype for relation-kind dispatch tokens.

Concrete subtypes — [`Membership`](@ref), [`Proximity`](@ref), [`KNN`](@ref),
[`Expression`](@ref) — are passed to `analyze` to select the algorithm, and
stored in the resulting `SpatialRelation` to enable re-dispatch.
"""
abstract type RelationKind end

"""
    Membership(; strict=false)
    Membership{strict}

Assignment of source observations to containing destination shapes.

- `strict=false` (default): containment is tested by point-in-polygon for
  transcripts, or centroid-in-polygon for cell shapes.
- `strict=true`: full geometric containment is required.

# See also
[`analyze`](@ref), [`SpatialRelation`](@ref)
"""
struct Membership{strict} <: RelationKind end
Membership(; strict::Bool=false) = Membership{strict}()

"""
    Proximity()

Relation kind for a pairwise distance matrix between elements.

The `weights` field of the resulting `SpatialRelation` is a symmetric N×N
`Float32` matrix of centroid distances. Boolean within-radius queries can be
derived as `weights .< r`.

# See also
[`distances`](@ref), [`KNN`](@ref)
"""
struct Proximity <: RelationKind end

"""
    KNN(; k=30)

Relation kind for a k-nearest-neighbour graph.

The resulting `SpatialRelation` has a N×k weight matrix and a flat N*k
`dst_ids` edge list. Requires a NearestNeighbors.jl backend:
`using NearestNeighbors`.

# See also
[`Proximity`](@ref), [`analyze`](@ref)
"""
struct KNN <: RelationKind
    k :: Int
end
KNN(; k::Int = 30) = KNN(k)

"""
    Expression()

Relation kind for a cell × gene transcript count matrix.

Produced by `analyze(Expression(), pts, cells)`, which assigns each transcript
to the cell containing it and accumulates counts. The `weights` field is a
`Float32` matrix of shape (n_cells × n_genes); `var` holds a named tuple with
a `:name` column of gene names.

# See also
[`analyze`](@ref), [`nobs`](@ref), [`nvar`](@ref), [`var_names`](@ref)
"""
struct Expression <: RelationKind end

# ── SpatialRelation ───────────────────────────────────────────────────────────

"""
    SpatialRelation{K<:RelationKind, W}

Weighted relation between two named spatial elements.

The relation kind `K` determines the semantics: `Expression` is a bipartite
cell × gene count matrix; `Membership` is a source-to-destination assignment;
`Proximity` and `KNN` are graph structures.

- `src`, `dst`: element names in the parent dataset
- `src_ids`, `dst_ids`: `instance_id` vectors identifying the rows/nodes
- `weights`: the relation data (`Matrix{Float32}` or `nothing`)
- `obs`: per-row metadata (Tables.jl-compatible)
- `var`: per-column metadata (for `Expression`: gene names via `:name`)
- `kind`: the `RelationKind` singleton

# Constructors
    SpatialRelation(kind, src, dst, src_ids, dst_ids, weights=nothing; obs)
    SpatialRelation(Expression(), src, src_ids, weights; obs, var)

# See also
[`analyze`](@ref), [`annotate`](@ref), [`nobs`](@ref), [`nvar`](@ref), [`var_names`](@ref)
"""
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

"""
    nobs(rel) → Int

Return the number of source observations (rows) in a `SpatialRelation`.
"""
nobs(rel::SpatialRelation) = length(rel.src_ids)

"""
    nvar(rel) → Int

Return the number of variables (columns) in an `Expression` relation.
Returns `0` for other relation kinds.

# See also
[`var_names`](@ref), [`nobs`](@ref)
"""
function nvar(rel::SpatialRelation{Expression})
    isnothing(rel.weights) ? 0 : size(rel.weights, 2)
end

"""
    var_names(rel) → Vector{String}

Return gene or variable names for an `Expression` relation.

Uses the `:name` column from `rel.var` if present; otherwise returns
string-formatted column indices.

# See also
[`nvar`](@ref), [`annotate`](@ref)
"""
function var_names(rel::SpatialRelation{Expression})
    t = Tables.columns(rel.var)
    Tables.columnnames(t) !== () && hasproperty(t, :name) ?
        collect(t.name) : string.(1:nvar(rel))
end

# ── annotate — pure: returns new SpatialRelation with obs column added ────────

"""
    annotate(rel, labels; key::Symbol) → SpatialRelation

Return a new `SpatialRelation` with `labels` added as column `key` in `rel.obs`.

Pure — does not modify `rel`. `labels` must have length equal to `nobs(rel)`.

```julia
types = assign_cell_types(rel)
rel2  = annotate(rel, types; key=:cell_type)
```

# See also
[`SpatialRelation`](@ref), [`nobs`](@ref)
"""
function annotate(rel::SpatialRelation, labels::AbstractVector; key::Symbol)
    length(labels) == nobs(rel) ||
        error("labels length ($(length(labels))) ≠ nobs(rel) ($(nobs(rel)))")
    new_obs = merge(Tables.columntable(rel.obs), NamedTuple{(key,)}((labels,)))
    SpatialRelation(rel.src, rel.dst, rel.src_ids, rel.dst_ids,
                    rel.weights, new_obs, rel.var, rel.kind)
end
