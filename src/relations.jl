# ── RelationKind hierarchy ────────────────────────────────────────────────────

"""
    RelationKind

Abstract supertype for relation-kind dispatch tokens.

Concrete subtypes — [`Membership`](@ref) and [`Expression`](@ref) — are passed
to `analyze` to select the algorithm, and
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
cell × gene count matrix; `Membership` is a source-to-destination assignment.

- `src`, `dst`: element names in the parent dataset
- `src_ids`, `dst_ids`: identifiers for the related observations. Point
  membership stores one-based point-row positions in `src_ids`; shape
  membership and expression relations use shape `instance_id` values.
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

# Non-expression relations have no variable metadata.
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
    source_ids(rel) -> Vector{Int32}

Return the source observation IDs stored by a relation.

For point [`Membership`](@ref), these are one-based row positions in the source
`SpatialPoints`. For shape membership and [`Expression`](@ref), they are source
shape `instance_id` values.
"""
source_ids(rel::SpatialRelation) = rel.src_ids

"""
    destination_ids(rel) -> Vector{Int32}

Return the destination IDs stored by a relation. Expression relations return
an empty vector because their columns are variables rather than destination
spatial objects.
"""
destination_ids(rel::SpatialRelation) = rel.dst_ids

"""
    nobs(rel) → Int

Return the number of relation rows. For `Expression`, this is the number of
source observations. For `Membership`, it is the number of matched pairs and
can exceed the number of unique sources when destination shapes overlap.
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
    obs_names(rel) → Vector{String}

Return observation names for a relation.

Uses the `:name` column from `rel.obs` if available; otherwise falls back to
string-formatted `src_ids`. Add names via `annotate(rel, names; key=:name)`.

# See also
[`var_names`](@ref), [`annotate`](@ref)
"""
function obs_names(rel::SpatialRelation)
    t = Tables.columntable(rel.obs)
    hasproperty(t, :name) ? collect(String, t.name) : string.(rel.src_ids)
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

# ── getindex — Expression row/col access by instance_id / gene name ──────────

const _RowIdx = Union{Colon, Integer, AbstractString,
                       AbstractVector{<:Integer}, AbstractVector{<:AbstractString}}
const _ColIdx = Union{Colon, AbstractString, AbstractVector{<:AbstractString}}

_row_idx(::SpatialRelation, ::Colon) = Colon()
function _row_idx(rel::SpatialRelation, id::Integer)
    i = findfirst(==(Int32(id)), rel.src_ids)
    isnothing(i) && error("instance_id $id not found in relation")
    i
end
function _row_idx(rel::SpatialRelation, name::AbstractString)
    obs = Tables.columntable(rel.obs)
    hasproperty(obs, :name) ||
        error("relation has no obs.name column; annotate with key=:name first")
    i = findfirst(==(name), obs.name)
    isnothing(i) && error("obs name $(repr(name)) not found; available: $(obs.name)")
    i
end
_row_idx(rel::SpatialRelation, ids::AbstractVector{<:Integer})     = [_row_idx(rel, id)   for id   in ids]
_row_idx(rel::SpatialRelation, names::AbstractVector{<:AbstractString}) = [_row_idx(rel, n) for n in names]

_col_idx(::SpatialRelation, ::Colon) = Colon()
function _col_idx(rel::SpatialRelation{Expression}, name::AbstractString)
    i = findfirst(==(name), var_names(rel))
    isnothing(i) && error("variable $(repr(name)) not found; available: $(var_names(rel))")
    i
end
_col_idx(rel::SpatialRelation{Expression}, names::AbstractVector{<:AbstractString}) =
    [_col_idx(rel, n) for n in names]

"""
    rel[rows, cols]

Index into an `Expression` relation. Dim 1 is the **container** (observation),
dim 2 is the **variable** (gene, cell type, etc.). Position determines semantics —
there is no ambiguity even when both axes carry string labels.

Both axes accept a singleton, a vector, or `:` (all):
- `rel[:, "Col1a1"]`               → `Vector{Float32}` — all obs, one gene
- `rel[47, :]`                     → `Vector{Float32}` — one obs, all genes
- `rel[47, "Col1a1"]`              → `Float32` scalar
- `rel[[47, 83], :]`               → `Matrix{Float32}`
- `rel[:, ["Col1a1","Epcam"]]`     → `Matrix{Float32}`

**Row axis** accepts integer instance IDs or, after `annotate(rel, names; key=:name)`,
string observation names:
- `rel["germinal_center", "Col1a1"]` — looks up row by `obs.name`
- `rel[["GC","TZ"], :]`             — multi-name row slice

**Column axis** always uses gene/variable names (strings from `var_names(rel)`).

# See also
[`obs_names`](@ref), [`var_names`](@ref), [`annotate`](@ref), [`SpatialRelation`](@ref)
"""
Base.getindex(rel::SpatialRelation{Expression}, rows::_RowIdx, cols::_ColIdx) =
    rel.weights[_row_idx(rel, rows), _col_idx(rel, cols)]

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
