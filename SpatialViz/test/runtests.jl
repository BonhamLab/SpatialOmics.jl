using Test
using CairoMakie
using SpatialViz
using SpatialOmicsBase
using DataFrames
using GeometryBasics: Polygon, Point2f, Circle

# ---------------------------------------------------------------------------
# Minimal fixtures — no I/O required
# ---------------------------------------------------------------------------

function _make_image(; nc=2, ny=32, nx=48)
    data = rand(Float32, nc, ny, nx)
    ax   = (c=nothing, y=nothing, x=nothing)
    meta = Dict{String,Any}("channel_labels" => ["DAPI", "PolyT"][1:nc])
    SpatialImage(data, ax, meta)
end

function _make_points(; n=20)
    coords   = rand(Float32, n, 2) .* 100f0
    features = DataFrame(gene = ["g$i" for i in 1:n])
    SpatialPoints(coords, features, Dict{String,Any}())
end

function _make_shapes()
    # Closed rings — first point repeated as last; required by GeometryOps centroid.
    polys = Vector{Any}([
        Polygon([Point2f(i, j), Point2f(i+2, j), Point2f(i+1, j+2), Point2f(i, j)])
        for i in 0:5:15 for j in 0:5:15
    ])
    n = length(polys)
    SpatialShapes(polys, DataFrame(id=1:n), Dict{String,Any}())
end

function _make_labels(; ny=32, nx=48)
    data = rand(0:5, ny, nx)
    SpatialLabels(data, Dict{String,Any}())
end

function _make_dataset()
    ds = spatial_dataset()
    ds.images["img"]   = _make_image()
    ds.points["pts"]   = _make_points()
    ds.shapes["cells"] = _make_shapes()
    ds.labels["lbl"]   = _make_labels()
    ds
end

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@testset "SpatialViz.jl" begin

    @testset "Exported symbols" begin
        for sym in [:SpatialExtent, :SpatialElementView, :SpatialDatasetView,
                    :extent, :intersects, :ImagePyramidSampler,
                    :spatial_panel, :composite]
            @test isdefined(SpatialViz, sym)
        end
    end

    @testset "heatmap convert_arguments — SpatialImage" begin
        img  = _make_image()
        args = Makie.convert_arguments(Makie.Heatmap, img)
        @test args isa Tuple
        @test length(args) == 3
    end

    @testset "heatmap convert_arguments — SpatialLabels" begin
        lbl  = _make_labels()
        args = Makie.convert_arguments(Makie.Heatmap, lbl)
        @test args isa Tuple
        @test length(args) == 3
    end

    @testset "scatter convert_arguments — SpatialPoints" begin
        pts      = _make_points()
        args     = Makie.convert_arguments(Makie.Scatter, pts)
        @test args isa Tuple
        pts_out  = first(args)
        @test pts_out isa AbstractVector
        @test length(pts_out) == 20
    end

    @testset "poly convert_arguments — SpatialShapes" begin
        shp  = _make_shapes()
        args = Makie.convert_arguments(Makie.Poly, shp)
        @test args isa Tuple
    end

    @testset "composite — no pyramid, fluor mode" begin
        img = _make_image(; nc=2)
        rgb = composite(img)
        @test eltype(rgb) == Makie.RGBf
        ny, nx = size(img.data, 2), size(img.data, 3)
        @test size(rgb) == (ny, nx)
    end

    @testset "composite — no pyramid, RGB mode" begin
        data = rand(Float32, 3, 16, 24)
        ax   = (c=nothing, y=nothing, x=nothing)
        meta = Dict{String,Any}("channel_labels" => ["R", "G", "B"])
        img  = SpatialImage(data, ax, meta)
        rgb  = composite(img)
        @test eltype(rgb) == Makie.RGBf
        @test size(rgb) == (16, 24)
    end

    @testset "composite — explicit channel selection" begin
        img = _make_image(; nc=2)
        rgb = composite(img; channel=[1], color=[:red])
        @test eltype(rgb) == Makie.RGBf
    end

    @testset "ImagePyramidSampler — size is positive" begin
        img     = _make_image(; nc=1, ny=32, nx=48)
        s       = ImagePyramidSampler(img, 1)
        nx_s, ny_s = size(s)
        @test nx_s > 0
        @test ny_s > 0
    end

    @testset "SpatialExtent from Rect2" begin
        r = Makie.Rect2f(10f0, 20f0, 30f0, 40f0)   # origin (10,20), widths (30,40)
        e = SpatialExtent(r)
        @test e.xmin ≈ 10.0
        @test e.xmax ≈ 40.0   # 10 + 30
        @test e.ymin ≈ 20.0
        @test e.ymax ≈ 60.0   # 20 + 40
    end

    @testset "view(ds, extent) returns SpatialDatasetView" begin
        ds  = _make_dataset()
        ext = SpatialExtent(0.0, 50.0, 0.0, 50.0)
        v   = view(ds, ext)
        @test v isa SpatialDatasetView
    end

    @testset "view(el, extent) returns SpatialElementView" begin
        pts = _make_points()
        ext = SpatialExtent(0.0, 50.0, 0.0, 50.0)
        v   = view(pts, ext)
        @test v isa SpatialElementView
    end

    @testset "scatter convert_arguments — SpatialElementView{SpatialPoints}" begin
        pts  = _make_points()
        ext  = SpatialExtent(0.0, 100.0, 0.0, 100.0)
        v    = view(pts, ext)
        args = Makie.convert_arguments(Makie.Scatter, v)
        @test args isa Tuple
    end

    @testset "poly convert_arguments — SpatialElementView{SpatialShapes}" begin
        shp  = _make_shapes()
        ext  = SpatialExtent(0.0, 20.0, 0.0, 20.0)
        v    = view(shp, ext)
        args = Makie.convert_arguments(Makie.Poly, v)
        @test args isa Tuple
    end

    @testset "spatial_panel — renders without error" begin
        ds = _make_dataset()
        fig, ax = spatial_panel(ds;
                                image_key="img", channel=1,
                                points_key="pts",
                                shapes_key="cells")
        @test fig isa Figure
        @test ax isa Axis
    end

    @testset "spatial_panel — with view" begin
        ds  = _make_dataset()
        ext = SpatialExtent(0.0, 30.0, 0.0, 30.0)
        v   = view(ds, ext)
        fig, ax = spatial_panel(ds;
                                image_key="img", channel=1,
                                points_key="pts",
                                view=v)
        @test fig isa Figure
    end

    @testset "heatmap! — SpatialImage" begin
        img = _make_image()
        fig = Figure()
        ax  = Axis(fig[1,1])
        @test_nowarn heatmap!(ax, img; channel=1, colormap=:grays)
    end

    @testset "heatmap! — SpatialLabels" begin
        lbl = _make_labels()
        fig = Figure()
        ax  = Axis(fig[1,1])
        @test_nowarn heatmap!(ax, lbl; colormap=:tab20)
    end

    @testset "scatter! — SpatialPoints" begin
        pts = _make_points()
        fig = Figure()
        ax  = Axis(fig[1,1])
        @test_nowarn scatter!(ax, pts; markersize=2)
    end

    @testset "poly! — SpatialShapes" begin
        shp = _make_shapes()
        fig = Figure()
        ax  = Axis(fig[1,1])
        @test_nowarn poly!(ax, shp; color=:transparent, strokecolor=:white)
    end

end
