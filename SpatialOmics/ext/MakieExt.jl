module MakieExt

using Makie
using SpatialOmics
using Colors: Gray, Colorant
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

function _select_level(v::SpatialImageColorView{C}; max_dim::Int=2048) where C
    raw = if isempty(v.pyramid)
        v.data
    else
        result = v.pyramid[end]
        for lvl in reverse(v.pyramid)
            maximum(size(lvl)) <= max_dim && (result = lvl; break)
        end
        result
    end
    dense = Array(raw)                          # step 1: bulk zarr read
    while maximum(size(dense)) > max_dim
        dense = restrict(dense)                 # step 2: downsample in-memory
    end
    eltype(dense) <: Colorant && return dense   # pre-colored composite: done
    display = v.transform !== nothing ? v.transform.(dense) : dense   # step 3
    colorview(v.colorant, display)              # step 4
end

function _pixel_extent(v::SpatialImageColorView)
    xi  = something(findfirst(==(:x), v.axes), 1)
    yi  = something(findfirst(==(:y), v.axes), 2)
    nx  = size(v.data, xi)
    ny  = size(v.data, yi)
    o   = apply(v.pixel_to_cs, SVector(0.0, 0.0))
    c   = apply(v.pixel_to_cs, SVector(Float64(nx), Float64(ny)))
    Float64[o[1], c[1]], Float64[o[2], c[2]]
end

function Makie.convert_arguments(P::Type{<:Image}, v::SpatialImageColorView)
    disp             = _select_level(v)
    x_range, y_range = _pixel_extent(v)
    return (x_range, y_range, disp)
end

Makie.convert_arguments(P::Type{<:Image}, img::SpatialImage) =
    convert_arguments(P, colorview(Gray, img))

# ── SpatialLabels → Heatmap ───────────────────────────────────────────────────

Makie.convert_arguments(P::Type{<:Heatmap}, lbl::SpatialLabels) =
    convert_arguments(P, Array(lbl.data))

end
