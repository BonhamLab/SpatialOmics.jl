# Coordinate systems

## The problem: multiple pixel spaces

Instruments like CosMx and MERFISH acquire data field-of-view by field-of-view.
Each FOV has its own pixel coordinate space. A single slide may contain hundreds
of FOVs, each offset and (for some instruments) rotated relative to the others.
Stitching them into a global coordinate space requires a per-FOV affine
transformation.

A naïve approach would stitch immediately: transform all coordinates into a
single global frame at load time, then discard the FOV structure. This works for
small datasets but fails at scale — you cannot analyse a single FOV independently
once its coordinates have been collapsed into global space, and global-pixel
images are far too large to hold in memory.

## Named coordinate systems as a graph

SpatialOmics models coordinate spaces explicitly as a directed acyclic graph
(DAG). Each node is a `CoordinateSystem` with a name, axis labels, and units.
Each edge is an `AbstractTransformation` carrying `src` and `dst` coordinate
system names.

```
  "fov_1_pixel"  ──Affine──►  "fov_1"  ──Affine──►  "global"
  "fov_2_pixel"  ──Affine──►  "fov_2"  ──Affine──►  "global"
  ...
```

Every element (`SpatialPoints`, `SpatialShapes`, etc.) stores a `coord_system`
name. That name places the element at a node in the graph. Elements at different
nodes are not directly comparable — the package raises an error rather than
silently misaligning data.

## Transformation primitives

Affine transformations in 2-D are represented as 3×3 augmented matrices in
homogeneous coordinates. This representation lets rotation, scaling, shear, and
translation be encoded uniformly, and lets sequential transforms be fused by
matrix multiplication.

The constructor helpers — `translation`, `scaling`, `rotation`, `flip_y` — each
produce an `Affine` with explicit `src` and `dst` names:

```julia
t = compose(
    scaling(0.325, 0.325, "pixel", "fov_1"),   # pixel size in µm
    translation(1024.0, 768.0, "fov_1", "global"),
)
push!(ds, CoordinateSystem("pixel"))
push!(ds, CoordinateSystem("fov_1"))
push!(ds, CoordinateSystem("global"))
push!(ds, t)
```

`compose` fuses two `Affine` transforms into one (matrix product), or wraps
mixed types in a `Sequence`. The `src`/`dst` chain must be consistent —
`a.dst == b.src` is enforced.

## Path resolution

`transform(ds, "pixel", "global")` calls `resolve`, which performs BFS over
the transform graph and returns a composed transformation from source to
destination. `Affine` edges are traversed in both directions (the inverse is
computed automatically); `Identity` edges are bidirectional by construction.

If the graph has no path between the requested systems, `resolve` raises an
error listing the known edges — a much clearer signal than a silent wrong
answer.

## Why not CoordinateTransformations.jl?

`CoordinateTransformations.jl` is a general-purpose library for function-based
transforms. SpatialOmics uses its own `Affine` type for two reasons: (1) the
augmented-matrix representation enables O(1) fusion via `compose`, which matters
when resolving paths through multi-hop graphs at load time; (2) every
transformation carries explicit `src` and `dst` names, making the graph
structure first-class rather than implicit in calling code.
