# make_fixtures.jl — developer script to create committable test fixtures
#
# Run from the SpatialOmics.jl repo root:
#   julia --project=docs/heavy test/make_fixtures.jl /path/to/xenium.zarr /path/to/visium.zarr
#
# Requires the full datasets downloaded from the SpatialData datasets page:
#   XENIUM_SRC  — Xenium FFPE Human Lung Cancer example (xenium_ex.zarr)
#   VISIUM_SRC  — Visium HD Mouse Small Intestine (visium_ex.zarr)
#
# The script writes overview figures, native fixtures, and ROI figures. To
# select a different region, inspect the overview assets, update the coordinate
# constants below, and run it again.

using SpatialOmics
using CairoMakie

length(ARGS) == 2 || error(
    "usage: julia --project=docs/heavy test/make_fixtures.jl XENIUM_ZARR VISIUM_ZARR",
)

const XENIUM_SRC = abspath(ARGS[1])
const VISIUM_SRC = abspath(ARGS[2])
const XENIUM_OUT = joinpath(@__DIR__, "data", "xenium_small.zarr")
const VISIUM_OUT = joinpath(@__DIR__, "data", "visium_small.zarr")
const ASSETS     = joinpath(@__DIR__, "..", "docs", "src", "assets")

mkpath(ASSETS)

# ══════════════════════════════════════════════════════════════════════════════
# XENIUM
# ══════════════════════════════════════════════════════════════════════════════

isdir(XENIUM_SRC) || error("Xenium source not found: $XENIUM_SRC")

@info "Loading Xenium dataset…"
xen = read(SpatialDataZarr(), XENIUM_SRC)
@info "Loaded" keys(elements(xen))

# ── Overview figure ───────────────────────────────────────────────────────────

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

# ── Fixture region ────────────────────────────────────────────────────────────
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

    sub["cell_labels"] = labels(roi, "cell_labels")

    rm(XENIUM_OUT; recursive=true, force=true)
    save!(sub; path=XENIUM_OUT)
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
    close(sub)
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
                            coord_system=coord_system(shp));
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
const VIS_CS   = coord_system(shapes(vis, VIS_SHAPES))

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
    save!(sub; path=VISIUM_OUT)
    @info "Visium fixture written → $VISIUM_OUT"

    fig2 = Figure(size=(600, 600))
    ax2  = Axis(fig2[1, 1]; aspect=DataAspect(), yreversed=true,
                title="Visium HD ROI — 16µm bins")
    image!(ax2, scaleminmax(channel(images(roi, VIS_IMAGE), 1)))
    poly!(ax2,   shapes(roi, VIS_SHAPES); color=(:steelblue, 0.4), strokewidth=0.3)
    tightlimits!(ax2)
    save(joinpath(ASSETS, "visium_roi.png"), fig2)
    @info "Saved ROI figure → $(joinpath(ASSETS, "visium_roi.png"))"
    close(sub)
end

close(xen)
close(vis)
