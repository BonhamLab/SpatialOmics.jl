
# ── Coordinate systems ────────────────────────────────────────────────────────

"""
    CoordinateSystem(name; axes=(:x, :y), units=("µm", "µm"))

Named 2-D coordinate system — a node in the dataset's transform graph.

Every spatial element belongs to one coordinate system by name. The graph of
`CoordinateSystem` nodes connected by `AbstractTransformation` edges is stored
in the parent `SpatialDataset`. Use `resolve` to find a composed path between
any two named systems.

# See also
[`SpatialDataset`](@ref), [`resolve`](@ref), [`AbstractTransformation`](@ref)
"""
struct CoordinateSystem
    name  :: String
    axes  :: NTuple{2, Symbol}
    units :: NTuple{2, String}
end

CoordinateSystem(name; axes=(:x, :y), units=("µm", "µm")) =
    CoordinateSystem(name, axes, units)

# ── Transformations ───────────────────────────────────────────────────────────

"""
    AbstractTransformation

Abstract supertype for all spatial coordinate transformations.

Concrete subtypes: [`Identity`](@ref), [`Affine`](@ref), [`Sequence`](@ref).
Every transformation carries `src` and `dst` coordinate system names.

# See also
[`apply`](@ref), [`apply!`](@ref), [`compose`](@ref), [`resolve`](@ref)
"""
abstract type AbstractTransformation end

"""
    Identity(src, dst)

Trivial transformation that returns its input unchanged.

Produced by `resolve` when `src == dst`, and used as the default `pixel_to_cs`
when no explicit transform is provided to `SpatialImage`.
"""
struct Identity <: AbstractTransformation
    src :: String
    dst :: String
end

"""
    Affine(matrix, src, dst)

2-D affine transformation stored as a 3×3 augmented matrix in homogeneous coordinates.

The matrix encodes rotation, scaling, shear, and translation in a single
`[R t; 0 0 1]` form, allowing sequential transforms to be fused by matrix
multiplication. Construct via [`translation`](@ref), [`scaling`](@ref),
[`rotation`](@ref), or [`flip_y`](@ref).

# See also
[`compose`](@ref), [`apply`](@ref)
"""
struct Affine <: AbstractTransformation
    matrix :: SMatrix{3, 3, Float64}   # augmented 2D: [R t; 0 0 1]
    src    :: String
    dst    :: String
end

"""
    Sequence(steps, src, dst)

Ordered composition of transformations applied left-to-right.

Produced by `resolve` when the path through the transform graph passes through
multiple intermediate coordinate systems, or when the steps cannot be fused
into a single `Affine` (e.g., if the path includes non-Affine steps).

# See also
[`compose`](@ref), [`resolve`](@ref)
"""
struct Sequence <: AbstractTransformation
    steps :: Vector{AbstractTransformation}
    src   :: String
    dst   :: String
end

# ── Constructor helpers ───────────────────────────────────────────────────────

"""
    translation(tx, ty, src, dst) → Affine

Affine transformation that shifts coordinates by `(tx, ty)`.
"""
function translation(tx::Real, ty::Real, src::String, dst::String)
    # SMatrix fills column-by-column; layout below reads left-to-right per row
    Affine(SMatrix{3,3,Float64}(1, 0, 0,
                                0, 1, 0,
                                tx, ty, 1), src, dst)
end

"""
    scaling(sx, sy, src, dst) → Affine

Affine transformation that scales the x-axis by `sx` and y-axis by `sy`.
"""
function scaling(sx::Real, sy::Real, src::String, dst::String)
    Affine(SMatrix{3,3,Float64}(sx, 0, 0,
                                0, sy, 0,
                                0,  0, 1), src, dst)
end

"""
    rotation(θ, src, dst) → Affine

Affine transformation for counter-clockwise rotation by angle `θ` (radians).
"""
function rotation(θ::Real, src::String, dst::String)
    c, s = cos(θ), sin(θ)
    # CCW rotation; SMatrix fills column-by-column, so cols are [c,s,0], [-s,c,0], [0,0,1]
    # giving matrix rows [c,-s,0], [s,c,0], [0,0,1]
    Affine(SMatrix{3,3,Float64}(c, s, 0,
                               -s, c, 0,
                                0, 0, 1), src, dst)
end

"""
    flip_y(src, dst) → Affine

Affine transformation that negates the y-axis (reflects about the x-axis).

Used to convert between image pixel space (y increases downward) and physical
space (y increases upward).
"""
function flip_y(src::String, dst::String)
    Affine(SMatrix{3,3,Float64}(1,  0, 0,
                                0, -1, 0,
                                0,  0, 1), src, dst)
end

"""
    compose(a, b, ...) → Affine or Sequence

Compose two or more transformations into a single transformation applied left-to-right.

When all arguments are `Affine`, the result is a fused `Affine` (matrix product).
Mixed types produce a `Sequence`. Raises an error if adjacent `src`/`dst` names
do not chain (`a.dst ≠ b.src`).

```julia
t = compose(scaling(0.325, 0.325, "pixel", "fov"), translation(1000.0, 500.0, "fov", "global"))
```

# See also
[`resolve`](@ref), [`apply`](@ref)
"""
function compose(a::Affine, b::Affine)
    a.dst == b.src || error("Cannot compose: $(a.dst) → $(b.src) mismatch")
    Affine(b.matrix * a.matrix, a.src, b.dst)
end

# ── Apply ─────────────────────────────────────────────────────────────────────

"""
    apply(t, data) → same type as data

Apply transformation `t` to `data`, returning a transformed copy.

Accepts a single `Point` / `SVector{2}`, a `Vector` of points, an N×2
`Matrix`, `SpatialPoints`, or `SpatialShapes`. The coordinate system name is
updated to `t.dst`.

# See also
[`apply!`](@ref), [`compose`](@ref), [`resolve`](@ref)
"""
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

"""
    resolve(transforms, src, dst) → AbstractTransformation

Find a transformation path from coordinate system `src` to `dst` through the
transform graph, and return the composed result.

Performs BFS over the directed graph of `(t.src → t.dst)` edges, treating
`Affine` edges as bidirectional (the inverse is computed automatically). When
all steps are `Affine`, they are fused into a single `Affine`; otherwise a
`Sequence` is returned. Raises an error if no path exists.

The two-argument convenience form `transform(ds, src, dst)` calls this function
using the dataset's transform list.

# See also
[`transform`](@ref), [`compose`](@ref), [`Sequence`](@ref)
"""
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
