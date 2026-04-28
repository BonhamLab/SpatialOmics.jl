module MakieExt

using Makie
using SpatialOmics

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

# ── SpatialImage → Heatmap ────────────────────────────────────────────────────
# Axes are stored as (x, y) for 2D or (x, y, c) for 3D after Zarr.jl's
# C→F reversal. Heatmap expects matrix[i,j] plotted at (i, j) — so data
# already needs no reorder for (x, y) slices; we only need Array() to
# materialise lazy zarr arrays.
#
# Default dispatches on channel 1 for 3D images.  Users wanting a specific
# channel can pass `img.data[:, :, ch]` directly to heatmap.

function Makie.convert_arguments(P::Type{<:Heatmap}, img::SpatialImage)
    N = ndims(img.data)
    if N == 2
        return convert_arguments(P, Array(img.data))
    else
        # find channel axis index
        ci = findfirst(==(:c), img.axes)
        ci === nothing && error("SpatialImage has no :c axis")
        slices = ntuple(d -> d == ci ? 1 : Colon(), N)
        return convert_arguments(P, Array(img.data[slices...]))
    end
end

# Coarsest pyramid level — fast overview without loading the full array
function _coarsest(img::SpatialImage)
    isempty(img.pyramid) ? img.data : img.pyramid[end]
end

# ── SpatialLabels → Heatmap ───────────────────────────────────────────────────

Makie.convert_arguments(P::Type{<:Heatmap}, lbl::SpatialLabels) =
    convert_arguments(P, Array(lbl.data))

end
