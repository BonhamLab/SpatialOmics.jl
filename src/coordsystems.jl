
# ── Coordinate systems ────────────────────────────────────────────────────────

struct CoordinateSystem
    name  :: String
    axes  :: NTuple{2, Symbol}
    units :: NTuple{2, String}
end

CoordinateSystem(name; axes=(:x, :y), units=("µm", "µm")) =
    CoordinateSystem(name, axes, units)

# ── Transformations ───────────────────────────────────────────────────────────

abstract type AbstractTransformation end

struct Identity <: AbstractTransformation
    src :: String
    dst :: String
end

struct Affine <: AbstractTransformation
    matrix :: SMatrix{3, 3, Float64}   # augmented 2D: [R t; 0 0 1]
    src    :: String
    dst    :: String
end

struct Sequence <: AbstractTransformation
    steps :: Vector{AbstractTransformation}
    src   :: String
    dst   :: String
end

# ── Constructor helpers ───────────────────────────────────────────────────────

function translation(tx::Real, ty::Real, src::String, dst::String)
    # SMatrix fills column-by-column; layout below reads left-to-right per row
    Affine(SMatrix{3,3,Float64}(1, 0, 0,
                                0, 1, 0,
                                tx, ty, 1), src, dst)
end

function scaling(sx::Real, sy::Real, src::String, dst::String)
    Affine(SMatrix{3,3,Float64}(sx, 0, 0,
                                0, sy, 0,
                                0,  0, 1), src, dst)
end

function rotation(θ::Real, src::String, dst::String)
    c, s = cos(θ), sin(θ)
    # CCW rotation; SMatrix fills column-by-column, so cols are [c,s,0], [-s,c,0], [0,0,1]
    # giving matrix rows [c,-s,0], [s,c,0], [0,0,1]
    Affine(SMatrix{3,3,Float64}(c, s, 0,
                               -s, c, 0,
                                0, 0, 1), src, dst)
end

function flip_y(src::String, dst::String)
    Affine(SMatrix{3,3,Float64}(1,  0, 0,
                                0, -1, 0,
                                0,  0, 1), src, dst)
end

function compose(a::Affine, b::Affine)
    a.dst == b.src || error("Cannot compose: $(a.dst) → $(b.src) mismatch")
    Affine(b.matrix * a.matrix, a.src, b.dst)
end

# ── Apply ─────────────────────────────────────────────────────────────────────

# Identity — return input unchanged for any supported container
apply(::Identity, pts::AbstractMatrix) = pts
apply(::Identity, p::StaticVector{2})  = p
apply(::Identity, pts::AbstractVector{<:StaticVector{2}}) = pts

# Affine on a single 2D point (covers GeometryBasics.Point2f, SVector{2}, etc.)
function apply(t::Affine, p::StaticVector{2,T}) where T<:Real
    aug = SVector{3,Float64}(p[1], p[2], 1.0)
    out = t.matrix * aug            # SMatrix{3,3,Float64} * SVector{3} → SVector{3}
    SVector(out[1], out[2])
end

# Affine on N×2 matrix (bulk, for internal array paths)
function apply(t::Affine, pts::AbstractMatrix{<:Real})
    n = size(pts, 1)
    aug = hcat(pts, ones(n))
    out = (t.matrix * aug')'
    out[:, 1:2]
end

# Vector of points — maps single-point apply across the collection
apply(t::AbstractTransformation, pts::AbstractVector{<:StaticVector{2}}) =
    map(p -> apply(t, p), pts)

# Sequence — generic fold; works for Matrix, Vector{SVector}, or any supported type
function apply(t::Sequence, data)
    foldl((d, step) -> apply(step, d), t.steps; init=data)
end

# Resolve ambiguity: Sequence × Vector{SVector} more specific than either above
apply(t::Sequence, pts::AbstractVector{<:StaticVector{2}}) =
    foldl((d, step) -> apply(step, d), t.steps; init=pts)

# ── Graph search to compose transforms on demand ──────────────────────────────

function _find_path(
    transforms::Vector{<:AbstractTransformation},
    src::String, dst::String
)
    src == dst && return [Identity(src, dst)]
    adj = Dict{String, Vector{Tuple{String, AbstractTransformation}}}()
    for t in transforms
        push!(get!(adj, t.src, []), (t.dst, t))
        # if Affine, also add inverse direction
        if t isa Affine
            push!(get!(adj, t.dst, []), (t.src, Affine(inv(Matrix(t.matrix)) |> SMatrix{3,3,Float64}, t.dst, t.src)))
        elseif t isa Identity
            push!(get!(adj, t.dst, []), (t.src, Identity(t.dst, t.src)))
        end
    end
    # BFS
    queue = [(src, AbstractTransformation[])]
    visited = Set{String}([src])
    while !isempty(queue)
        (cur, path) = popfirst!(queue)
        for (nxt, t) in get(adj, cur, [])
            nxt in visited && continue
            newpath = vcat(path, [t])
            nxt == dst && return newpath
            push!(visited, nxt)
            push!(queue, (nxt, newpath))
        end
    end
    nothing
end

function resolve(
    transforms::Vector{<:AbstractTransformation},
    src::String, dst::String
)
    path = _find_path(transforms, src, dst)
    path === nothing && error(
        "No transform path from \"$src\" to \"$dst\". Known transforms: " *
        join(["$(t.src)→$(t.dst)" for t in transforms], ", ")
    )
    length(path) == 1 && return path[1]
    # collapse Affines where possible
    if all(t -> t isa Affine, path)
        return foldl(compose, path)
    end
    Sequence(path, src, dst)
end
