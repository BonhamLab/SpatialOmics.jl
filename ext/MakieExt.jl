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

# ── SpatialImage/SpatialImageColorView → Image (zoom-responsive) ─────────────
# For images with pyramid levels, image!(ax, ...) pushes a new Observable value
# when the axis zoom changes, selecting the pyramid level whose full-res/screen-px
# ratio best matches the current viewport. Standard Makie reactive pattern.

function _materialise_level(v::SpatialImageColorView{C}, lvl) where C
    xi = something(findfirst(==(:x), v.axes), 1)
    yi = something(findfirst(==(:y), v.axes), 2)
    dense = Array(lvl)
    eltype(dense) <: Colorant && return yi < xi ? permutedims(dense, (2, 1)) : dense
    disp = v.transform !== nothing ? v.transform.(dense) : _to_colorable(dense)
    cv   = colorview(C, disp)
    yi < xi ? permutedims(collect(cv), (2, 1)) : collect(cv)
end

function Makie.image!(ax::Makie.Axis, img::SpatialImage; kw...)
    Makie.image!(ax, colorview(Gray, img); kw...)
end

function Makie.image!(ax::Makie.Axis, v::SpatialImageColorView; kw...)
    x_range, y_range = _pixel_extent(v)
    xi = something(findfirst(==(:x), v.axes), 1)
    yi = something(findfirst(==(:y), v.axes), 2)

    if isempty(v.pyramid)
        d   = _select_level(v)
        out = ndims(d) == 2 && yi < xi ? permutedims(d, (2, 1)) : d
        return Makie.image!(ax, x_range, y_range, out; kw...)
    end

    pyr_scales = [round(Int, size(v.data, yi) / size(l, yi)) for l in v.pyramid]

    cur_lev = Ref(lastindex(v.pyramid))
    img_obs = Observable(_materialise_level(v, v.pyramid[end]))
    plt     = Makie.image!(ax, x_range, y_range, img_obs; kw...)

    on(ax.finallimits; update=true) do lims
        ax.scene.viewport[].widths[1] == 0 && return
        vis_w  = max(lims.widths[1], 1.0)
        scn_w  = max(ax.scene.viewport[].widths[1], 1)
        target = max(1, floor(Int, vis_w / scn_w))

        best = 1
        for i in lastindex(pyr_scales):-1:1
            pyr_scales[i] <= target && (best = i; break)
        end

        best == cur_lev[] && return
        cur_lev[] = best
        img_obs[] = _materialise_level(v, v.pyramid[best])
    end

    plt
end

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

# ── Interactive ROI selection ─────────────────────────────────────────────────

# SpatialExtent from current axis limits
function SpatialOmics.SpatialExtent(ax::Makie.Axis; coord_system::String="")
    r = ax.finallimits[]
    o = minimum(r); w = widths(r)
    SpatialExtent(o[1], o[1]+w[1], o[2], o[2]+w[2]; coord_system)
end

# Per-axis session state for in-progress polygon drawing
const _active_roi_sessions = Dict{Makie.Axis, NamedTuple}()

function _teardown_roi!(ax)
    sess = get(_active_roi_sessions, ax, nothing)
    sess === nothing && return
    delete!(ax, sess.preview_lines)
    delete!(ax, sess.preview_dots)
    delete!(ax, sess.first_dot)
    Observables.off(sess.h_mouse)
    Observables.off(sess.h_key)
    delete!(_active_roi_sessions, ax)
end

# Interactive polygon drawing. Left-click adds vertices; click within snap_px of
# the first vertex (or press Enter) to close. Escape cancels.
# Returns an Observable — check obs[] after closing.
function SpatialOmics.select(ax::Makie.Axis, ::Type{SpatialROI};
                              coord_system::String="", snap_px::Real=10, priority::Int=2)
    _teardown_roi!(ax)

    result      = Observable{Union{Nothing, SpatialROI}}(nothing)
    vertices    = Point2f[]
    preview_pts = Observable(Point2f[])
    first_pt    = Observable(Point2f[])

    preview_lines = lines!(ax,  preview_pts; color=(:red, 0.7), linewidth=2)
    preview_dots  = scatter!(ax, preview_pts; color=:red,  markersize=8)
    first_dot     = scatter!(ax, first_pt;   color=:cyan, markersize=14, marker=:circle)

    function close_polygon!()
        ring = copy(vertices)
        first(ring) ≈ last(ring) || push!(ring, ring[1])
        result[] = SpatialROI(Polygon(ring); coord_system)
        _teardown_roi!(ax)
    end

    function update_preview!()
        preview_pts[] = length(vertices) >= 2 ?
            vcat(vertices, [vertices[1]]) : copy(vertices)
        first_pt[] = isempty(vertices) ? Point2f[] : [vertices[1]]
    end

    h_mouse = on(events(ax.scene).mousebutton, priority=priority) do event
        is_mouseinside(ax.scene) || return Consume(false)
        event.action == Mouse.press && event.button == Mouse.left || return Consume(false)
        if length(vertices) >= 3
            first_screen = Makie.project(ax.scene, vertices[1])
            norm(first_screen - Point2f(events(ax.scene).mouseposition[])) < snap_px &&
                (close_polygon!(); return Consume(false))
        end
        push!(vertices, mouseposition(ax))
        update_preview!()
        return Consume(false)
    end

    h_key = on(events(ax.scene).keyboardbutton, priority=priority) do event
        event.action == Keyboard.press || return Consume(false)
        if event.key == Keyboard.enter || event.key == Keyboard.kp_enter
            length(vertices) >= 3 && close_polygon!()
        elseif event.key == Keyboard.escape
            result[] = nothing
            _teardown_roi!(ax)
        end
        return Consume(false)
    end

    _active_roi_sessions[ax] = (; preview_lines, preview_dots, first_dot, h_mouse, h_key)
    return result
end

end
