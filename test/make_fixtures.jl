# make_fixtures.jl — developer script to create committable test fixtures
#
# Run from the SpatialOmics.jl repo root:
#   julia --project=. test/make_fixtures.jl
#
# Requires the full datasets downloaded from the SpatialData datasets page:
#   XENIUM_SRC  — Xenium Mouse Brain example (xenium_ex.zarr)
#   VISIUM_SRC  — Visium HD Mouse Small Intestine (visium_ex.zarr)
#
# Workflow:
#   1. Run the script once — it saves overview figures and exits.
#   2. Open the overview PNGs (written to docs/src/assets/), pick a region.
#   3. Fill in the coordinate constants below (STEP 2 blocks).
#   4. Re-run — fixture zarrs are written to test/data/.

using SpatialOmics
using CairoMakie
using StaticArrays

const XENIUM_SRC = "/home/kevin/Repos/stx_dev/test_data/experiments/xenium_ex.zarr"
const VISIUM_SRC = "/home/kevin/Repos/stx_dev/test_data/experiments/visium_ex.zarr"
const XENIUM_OUT = joinpath(@__DIR__, "data", "xenium_small.zarr")
const VISIUM_OUT = joinpath(@__DIR__, "data", "visium_small.zarr")
const ASSETS     = joinpath(@__DIR__, "..", "docs", "src", "assets")

mkpath(ASSETS)

# ── Helpers ───────────────────────────────────────────────────────────────────

# Maps a physical-space point back to pixel coordinates using the image transform.
function _to_pixel(t::Affine, x::Real, y::Real)
    m_inv = inv(Matrix(t.matrix))
    v = m_inv * [Float64(x), Float64(y), 1.0]
    (v[1], v[2])
end
_to_pixel(::Identity, x::Real, y::Real) = (Float64(x), Float64(y))

# Clamps a pixel range to valid array bounds (1-based, inclusive).
_px_clamp(lo, hi, n) = (clamp(floor(Int, min(lo, hi)) + 1, 1, n),
                         clamp(ceil(Int,  max(lo, hi)),     1, n))

# Adjust pixel_to_cs for a crop: new pixel [1,1] = old pixel [xlo, ylo].
function _shift_origin(t::Affine, xlo::Int, ylo::Int)
    dx, dy = Float64(xlo - 1), Float64(ylo - 1)
    S = SMatrix{3,3,Float64}(1, 0, 0, 0, 1, 0, dx, dy, 1)
    Affine(t.matrix * S, t.src, t.dst)
end
_shift_origin(t::Identity, ::Int, ::Int) = t

# Crop a SpatialLabels to a pixel rectangle and filter the instance_map.
<<<<<<< HEAD
# axes are (:x, :y), so dim1=x, dim2=y — index as [xlo:xhi, ylo:yhi].
function crop_labels(lbl::SpatialLabels, ylo::Int, yhi::Int, xlo::Int, xhi::Int)
    raw = Array(lbl.data[xlo:xhi, ylo:yhi])
=======
function crop_labels(lbl::SpatialLabels, ylo::Int, yhi::Int, xlo::Int, xhi::Int)
    raw = Array(lbl.data[ylo:yhi, xlo:xhi])
>>>>>>> fix
    present = Set(raw)
    imap = Dict(k => v for (k, v) in lbl.instance_map if k in present)
    p2cs = _shift_origin(lbl.pixel_to_cs, xlo, ylo)
    SpatialLabels(raw; axes=lbl.axes, instance_map=imap,
                  coord_system=lbl.coord_system, pixel_to_cs=p2cs)
end

# ══════════════════════════════════════════════════════════════════════════════
# XENIUM
# ══════════════════════════════════════════════════════════════════════════════

isdir(XENIUM_SRC) || error("Xenium source not found: $XENIUM_SRC")

@info "Loading Xenium dataset…"
xen = read(SpatialDataZarr(), XENIUM_SRC)
@info "Loaded" keys(elements(xen))

# ── STEP 1: Overview figure ───────────────────────────────────────────────────
# Run this block first, inspect xenium_overview.png, then fill in coordinates.

let
    fig = Figure(size=(900, 900))
    ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true,
               title="Xenium — transcript overview (200k subsampled)")
    tx  = points(xen, "transcripts")
    n   = length(coords(tx))
    scatter!(ax, subsample(tx, min(n, 200_000)); markersize=1.5, color=(:black, 0.2))
    tightlimits!(ax)
    path = joinpath(ASSETS, "xenium_overview.png")
    save(path, fig)
    @info "Saved overview → $path — inspect to pick XMIN/XMAX/YMIN/YMAX (coordinates are in µm)"
end

# ── STEP 2: Fill in these values after viewing the overview ───────────────────
# Choose a ~200µm × 200µm region with good transcript density and visible cells.
# coord_system must match points(xen,"transcripts").coord_system.

const XEN_XMIN = 30000.0
const XEN_XMAX = 30200.0
const XEN_YMIN =  7500.0
const XEN_YMAX =  7700.0
const XEN_CS   = "global"  # verify with coord_systems(xen)

if XEN_XMAX == XEN_XMIN
    @info "Xenium extent not set — fill in XEN_XMIN/XMAX/YMIN/YMAX and re-run"
else
    @info "Creating Xenium fixture…"
    ext = SpatialExtent(XEN_XMIN, XEN_XMAX, XEN_YMIN, XEN_YMAX; coord_system=XEN_CS)
    roi = view(xen, ext)

    sub = SpatialDataset()
    for (_, cs) in xen.coord_systems; push!(sub, cs); end
    for t in xen.transforms;         push!(sub, t);  end

    sub["transcripts"]     = collect(points(roi, "transcripts"))
    sub["cell_boundaries"] = collect(shapes(roi, "cell_boundaries"))

    # Image: view() returns a lazily cropped SpatialImage with adjusted pixel_to_cs
    sub["morphology_focus"] = images(roi, "morphology_focus")

    # Labels: manual crop (view on SpatialDatasetView returns labels as-is)
    img_ref = images(xen, "morphology_focus")
    lbl     = labels(xen, "cell_labels")
    ny, nx  = size(lbl.data, findfirst(==(:y), lbl.axes)),
              size(lbl.data, findfirst(==(:x), lbl.axes))
    px1 = _to_pixel(img_ref.pixel_to_cs, XEN_XMIN, XEN_YMIN)
    px2 = _to_pixel(img_ref.pixel_to_cs, XEN_XMAX, XEN_YMAX)
    xlo, xhi = _px_clamp(px1[1], px2[1], nx)
    ylo, yhi = _px_clamp(px1[2], px2[2], ny)
    sub["cell_labels"] = crop_labels(lbl, ylo, yhi, xlo, xhi)

    rm(XENIUM_OUT; recursive=true, force=true)
    write!(sub, XENIUM_OUT, SpatialDataZarr())
    @info "Xenium fixture written → $XENIUM_OUT" size=Base.format_bytes(
        sum(filesize(f) for (r,_,fs) in walkdir(XENIUM_OUT) for f in joinpath.(r,fs)))

    # ROI composite figure
    fig2 = Figure(size=(600, 600))
    ax2  = Axis(fig2[1, 1]; aspect=DataAspect(), yreversed=true,
                title="Xenium ROI — $(round(Int, XEN_XMAX-XEN_XMIN)) × $(round(Int, XEN_YMAX-XEN_YMIN)) µm")
    image!(ax2, scaleminmax(channel(images(roi, "morphology_focus"), 1)))
    poly!(ax2,  shapes(roi, "cell_boundaries");  color=:transparent,
                                                   strokecolor=:cyan, strokewidth=0.5)
    scatter!(ax2, points(roi, "transcripts");      markersize=1.5, color=(:red, 0.4))
    tightlimits!(ax2)
    save(joinpath(ASSETS, "xenium_roi.png"), fig2)
    @info "Saved ROI figure → $(joinpath(ASSETS, "xenium_roi.png"))"
end

# ══════════════════════════════════════════════════════════════════════════════
# VISIUM
# ══════════════════════════════════════════════════════════════════════════════

isdir(VISIUM_SRC) || error("Visium source not found: $VISIUM_SRC")

@info "Loading Visium dataset…"
vis = read(SpatialDataZarr(), VISIUM_SRC)
@info "Loaded" keys(elements(vis))

# Key names (long because they come from the Python SpatialData store)
const VIS_SHAPES = "Visium_HD_Mouse_Small_Intestine_square_016um"
const VIS_IMAGE  = "Visium_HD_Mouse_Small_Intestine_lowres_image"

let
    fig = Figure(size=(900, 900))
    ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true,
               title="Visium HD — 16µm bins overview")
    shp = shapes(vis, VIS_SHAPES)
    n   = length(geometries(shp))
    idx = n > 5_000 ? rand(1:n, 5_000) : 1:n
    poly!(ax, SpatialShapes(geometries(shp)[idx]; instance_id=instance_id(shp)[idx],
                            coord_system=shp.coord_system);
          color=:steelblue, strokewidth=0)
    tightlimits!(ax)
    path = joinpath(ASSETS, "visium_overview.png")
    save(path, fig)
    @info "Saved overview → $path — inspect this figure, then fill in VIS_XMIN/XMAX/YMIN/YMAX below"
end

const VIS_XMIN = 3000.0
const VIS_XMAX = 3500.0
const VIS_YMIN = 2000.0
const VIS_YMAX = 2500.0
const VIS_CS   = shapes(vis, VIS_SHAPES).coord_system  # auto-detected from element

if VIS_XMAX == VIS_XMIN
    @info "Visium extent not set — fill in VIS_XMIN/XMAX/YMIN/YMAX and re-run"
else
    @info "Creating Visium fixture…"
    ext = SpatialExtent(VIS_XMIN, VIS_XMAX, VIS_YMIN, VIS_YMAX; coord_system=VIS_CS)
    roi = view(vis, ext)

    sub = SpatialDataset()
    for (_, cs) in vis.coord_systems; push!(sub, cs); end
    for t in vis.transforms;         push!(sub, t);  end

    sub[VIS_SHAPES] = collect(shapes(roi, VIS_SHAPES))
    sub[VIS_IMAGE]  = images(roi, VIS_IMAGE)

    rm(VISIUM_OUT; recursive=true, force=true)
    write!(sub, VISIUM_OUT, SpatialDataZarr())
    @info "Visium fixture written → $VISIUM_OUT"

    fig2 = Figure(size=(600, 600))
    ax2  = Axis(fig2[1, 1]; aspect=DataAspect(), yreversed=true,
                title="Visium HD ROI — 16µm bins")
    image!(ax2, scaleminmax(channel(images(roi, VIS_IMAGE), 1)))
    poly!(ax2,   shapes(roi, VIS_SHAPES); color=(:steelblue, 0.4), strokewidth=0.3)
    tightlimits!(ax2)
    save(joinpath(ASSETS, "visium_roi.png"), fig2)
    @info "Saved ROI figure → $(joinpath(ASSETS, "visium_roi.png"))"
end
