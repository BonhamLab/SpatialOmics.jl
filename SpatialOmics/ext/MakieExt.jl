module MakieExt

using Makie
using SpatialOmics
using Colors: Gray, Colorant
using FixedPointNumbers: FixedPoint, Normed
using StaticArrays: SVector
using ImageBase: restrict

# ── SpatialPoints → Scatter ───────────────────────────────────────────────────

Makie.convert_arguments(P::Type{<:Scatter}, pts::SpatialPoints) =
    convert_arguments(P, pts.coords)

Makie.convert_arguments(P::Type{<:Scatter}, v::SpatialElementView{<:SpatialPoints}) =
    convert_arguments(P, coords(v))

# ── SpatialShapes → Poly ──────────────────────────────────────────────────────

Makie.convert_arguments(P::Type{<:Poly}, shps::SpatialShapes) =
    convert_arguments(P, geometries(shps))

Makie.convert_arguments(P::Type{<:Poly}, v::SpatialElementView{<:SpatialShapes}) =
    convert_arguments(P, geometries(v))

# ── SpatialImageColorView → Image ─────────────────────────────────────────────
# Selects the finest pyramid level whose longest dimension ≤ max_dim, then:
#   1. Materialises zarr via Array() — bulk read (fast)
#   2. Restricts further if still oversized
#   3. Applies display transform (f.(dense)) over in-memory array (fast)
#   4. Applies colorview to produce Colorant array for Makie
# This ordering is critical: wrapping zarr in any lazy transform before Array()
# defeats the chunk-based bulk-read path and causes catastrophic slowdown.

_to_colorable(arr::AbstractArray) = arr
_to_colorable(arr::AbstractArray{T}) where {T<:Unsigned} =
    reinterpret(Normed{T, 8 * sizeof(T)}, arr)

function _select_level(v::SpatialImageColorView{C}; max_dim::Int=4096) where C
    raw = if isempty(v.pyramid) || maximum(size(v.data)) <= max_dim
        v.data                            # full-res fits (common for crops)
    else
        result = v.pyramid[end]           # fallback: coarsest
        for lvl in v.pyramid              # finest → coarsest; take first that fits
            maximum(size(lvl)) <= max_dim && (result = lvl; break)
        end
        result
    end
    dense = Array(raw)                          # step 1: bulk zarr read
    while maximum(size(dense)) > max_dim
        dense = restrict(dense)                 # step 2: downsample in-memory
    end
    eltype(dense) <: Colorant && return dense   # pre-colored composite: done
    display = v.transform !== nothing ? v.transform.(dense) : _to_colorable(dense)  # step 3
    colorview(v.colorant, display)              # step 4
end

function _pixel_extent(v::SpatialImageColorView)
    xi  = something(findfirst(==(:x), v.axes), 1)
    yi  = something(findfirst(==(:y), v.axes), 2)
    nx  = size(v.data, xi)
    ny  = size(v.data, yi)
    o   = apply(v.pixel_to_cs, SVector(0.0, 0.0))
    c   = apply(v.pixel_to_cs, SVector(Float64(nx), Float64(ny)))
    (o[1], c[1]), (o[2], c[2])
end

function Makie.convert_arguments(P::Type{<:Image}, v::SpatialImageColorView)
    disp             = _select_level(v)
    x_range, y_range = _pixel_extent(v)
    # Makie image!(ax, x_range, y_range, data): data[i,j] → position (x[i], y[j]).
    # So data must be x-first (first dim = x, second = y).
    # If data is y-first (yi < xi), transpose to make it x-first.
    xi  = something(findfirst(==(:x), v.axes), 1)
    yi  = something(findfirst(==(:y), v.axes), 2)
    out = (ndims(disp) == 2 && yi < xi) ? permutedims(disp, (2, 1)) : disp
    return (x_range, y_range, out)
end

Makie.convert_arguments(P::Type{<:Image}, img::SpatialImage) =
    convert_arguments(P, colorview(Gray, img))

# ── SpatialLabels → Heatmap ───────────────────────────────────────────────────

Makie.convert_arguments(P::Type{<:Heatmap}, lbl::SpatialLabels) =
    convert_arguments(P, Array(lbl.data))

# ── PointDensity → Heatmap ────────────────────────────────────────────────────
# density(pts; resolution=512, feature="EPCAM") → PointDensity
# heatmap!(ax, density(pts; resolution=256)) bins transcripts into a 2D grid.

function Makie.convert_arguments(P::Type{<:Heatmap}, d::PointDensity)
    cds = d.feature !== nothing ? coords(d.pts, d.feature) : coords(d.pts)
    isempty(cds) && return (0f0:1f0, 0f0:1f0, zeros(Float32, 1, 1))
    xs = Float32[p[1] for p in cds]
    ys = Float32[p[2] for p in cds]
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    n  = d.resolution
    dx = max((xmax - xmin) / n, eps(Float32))
    dy = max((ymax - ymin) / n, eps(Float32))
    counts = zeros(Float32, n, n)
    for i in eachindex(xs)
        ix = clamp(ceil(Int, (xs[i] - xmin) / dx), 1, n)
        iy = clamp(ceil(Int, (ys[i] - ymin) / dy), 1, n)
        counts[ix, iy] += 1f0
    end
    (range(xmin, xmax; length=n), range(ymin, ymax; length=n), counts)
end

# ── ShapeColorView → Poly ─────────────────────────────────────────────────────
# poly!(ax, cells, rel; color_by=:cell_type, colormap=:tab20)
# Constructs a ShapeColorView and draws polygons colored by an obs column.

function Makie.poly!(ax::Makie.Axis, v::ShapeColorView; kw...)
    geoms = geometries(v.shapes)
    vals  = v.color_by ∈ propertynames(v.rel.obs) ?
            v.rel.obs[v.color_by] : ones(Int, length(geoms))
    uniq  = unique(vals)
    cmap  = Makie.to_colormap(v.colormap)
    colors = [cmap[mod1(findfirst(==(val), uniq), length(cmap))] for val in vals]
    Makie.poly!(ax, geoms; color=colors, kw...)
end

function Makie.poly!(ax::Makie.Axis, cells::SpatialShapes, rel::SpatialRelation;
                     color_by::Symbol=:label, colormap=:tab20, kw...)
    Makie.poly!(ax, ShapeColorView(cells, rel, color_by, colormap); kw...)
end

end
