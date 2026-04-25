using SpatialOmics
using Test

@testset "SpatialOmics M1" begin

    @testset "CoordinateSystem" begin
        cs = CoordinateSystem("global"; axes=(:x, :y), units=("µm", "µm"))
        @test cs.name == "global"
        @test cs.axes == (:x, :y)
        @test cs.units == ("µm", "µm")
    end

    @testset "Transformations — construction" begin
        t = translation(10.0, -5.0, "fov_1", "global")
        @test t isa Affine
        @test t.src == "fov_1"
        @test t.dst == "global"

        s = scaling(2.0, 2.0, "px", "µm")
        @test s isa Affine

        r = rotation(π/4, "a", "b")
        @test r isa Affine

        f = flip_y("local", "global")
        @test f isa Affine
    end

    @testset "Transformations — apply" begin
        # translation
        pts = [1.0 2.0; 3.0 4.0]   # 2×2
        t = translation(10.0, 20.0, "a", "b")
        out = apply(t, pts)
        @test out ≈ [11.0 22.0; 13.0 24.0]

        # flip_y
        f = flip_y("a", "b")
        out2 = apply(f, pts)
        @test out2 ≈ [1.0 -2.0; 3.0 -4.0]

        # identity
        id = Identity("a", "b")
        @test apply(id, pts) === pts

        # compose two translations
        t1 = translation(1.0, 0.0, "a", "b")
        t2 = translation(0.0, 1.0, "b", "c")
        tc = compose(t1, t2)
        @test apply(tc, [0.0 0.0]) ≈ [1.0 1.0]
    end

    @testset "Transformations — resolve_transform / Dijkstra" begin
        transforms = AbstractTransformation[
            translation(100.0, 200.0, "fov_1", "global"),
            scaling(0.5, 0.5, "px", "µm"),
        ]
        t = resolve_transform(transforms, "fov_1", "global")
        @test t isa Affine

        # no path
        @test_throws ErrorException resolve_transform(transforms, "nowhere", "global")

        # identity (same src == dst)
        t2 = resolve_transform(transforms, "global", "global")
        @test t2 isa Identity
    end

    @testset "Sequence apply" begin
        t1 = translation(1.0, 0.0, "a", "b")
        t2 = translation(0.0, 1.0, "b", "c")
        seq = Sequence([t1, t2], "a", "c")
        out = apply(seq, [0.0 0.0])
        @test out ≈ [1.0 1.0]
    end

    @testset "apply on SVector / Point2f" begin
        using StaticArrays
        t = translation(10.0, 20.0, "a", "b")

        # single SVector{2}
        p = SVector(1.0, 2.0)
        out = apply(t, p)
        @test out isa SVector{2}
        @test out ≈ SVector(11.0, 22.0)

        # Identity on SVector
        id = Identity("a", "b")
        @test apply(id, p) === p

        # vector of SVectors
        pts = [SVector(0.0, 0.0), SVector(1.0, 1.0), SVector(2.0, 3.0)]
        outs = apply(t, pts)
        @test outs isa Vector
        @test outs[1] ≈ SVector(10.0, 20.0)
        @test outs[3] ≈ SVector(12.0, 23.0)

        # Identity on vector of SVectors
        @test apply(id, pts) === pts

        # Sequence on vector of SVectors
        t1 = translation(1.0, 0.0, "a", "b")
        t2 = translation(0.0, 1.0, "b", "c")
        seq = Sequence([t1, t2], "a", "c")
        svec_pts = [SVector(0.0, 0.0)]
        @test apply(seq, svec_pts)[1] ≈ SVector(1.0, 1.0)

        # flip_y on SVector
        f = flip_y("a", "b")
        @test apply(f, SVector(3.0, 4.0)) ≈ SVector(3.0, -4.0)
    end

    @testset "BackingStore — tempdir" begin
        bs = BackingStore()
        @test isdir(bs.path)
        @test startswith(basename(bs.path), "spatialomics_")
        @test bs.owned == true
        zarr_json = joinpath(bs.path, "zarr.json")
        @test isfile(zarr_json)
        p = bs.path
        SpatialOmics._cleanup!(bs)
        @test !isdir(p)
    end

    @testset "BackingStore — user path" begin
        mktempdir() do d
            bs = BackingStore(; path=d)
            @test bs.owned == false
            @test isfile(joinpath(bs.path, "zarr.json"))
        end
    end

    @testset "SpatialDataset — construction + cleanup" begin
        ds = SpatialDataset()
        p = ds.backing.path
        @test isdir(p)
        @test ds.backing.owned == true
        close(ds)
        @test !isdir(p)
    end

    @testset "SpatialDataset — keep!" begin
        mktempdir() do d
            target = joinpath(d, "myds.zarr")
            ds = SpatialDataset()
            p = ds.backing.path
            keep!(ds, target)
            @test isdir(target)
            @test ds.backing.owned == false
            @test !isdir(p)         # scratch removed after copy
            close(ds)               # should be a no-op (owned=false)
            @test isdir(target)
        end
    end

    @testset "SpatialDataset — with_dataset" begin
        path_ref = Ref("")
        with_dataset() do ds
            path_ref[] = ds.backing.path
            @test isdir(ds.backing.path)
        end
        @test !isdir(path_ref[])
    end

    @testset "SpatialDataset — coord systems and transforms" begin
        ds = SpatialDataset()
        try
            push!(ds, CoordinateSystem("global"; units=("µm", "µm")))
            push!(ds, CoordinateSystem("fov_1"; axes=(:x, :y), units=("px", "px")))
            @test "global" in coord_systems(ds)
            @test "fov_1" in coord_systems(ds)

            t = translation(500.0, 300.0, "fov_1", "global")
            push!(ds, t)
            resolved = transform(ds, "fov_1", "global")
            @test resolved isa Affine

            pts = [0.0 0.0; 10.0 20.0]
            out = apply(resolved, pts)
            @test out ≈ [500.0 300.0; 510.0 320.0]
        finally
            close(ds)
        end
    end

    @testset "SpatialDataset — element setindex/getindex" begin
        ds = SpatialDataset()
        try
            ds["test"] = (x = 1, y = 2)
            @test haskey(ds, "test")
            @test ds["test"] == (x = 1, y = 2)
            @test "test" in collect(keys(ds))
        finally
            close(ds)
        end
    end

end

@testset "SpatialOmics M2" begin

    using GeometryBasics
    using GeoInterface

    @testset "SpatialPoints — bare constructor" begin
        pts = SpatialPoints([Point2f(1, 2), Point2f(3, 4), Point2f(5, 6)]; coord_system="global")
        @test pts isa SpatialPoints{Float32}
        @test length(pts) == 3
        @test coord_system(pts) == "global"
        @test coords(pts)[1] == Point2f(1, 2)
        @test features(pts) == String[]
        @test all(feature_ids(pts) .== 0)
        @test all(pts.instance_id .== 0)    # internal field
    end

    @testset "SpatialPoints — Tables constructor" begin
        pts = SpatialPoints(
            (x = [1.0f0, 2.0f0, 3.0f0],
             y = [4.0f0, 5.0f0, 6.0f0],
             gene = ["Actb", "Gapdh", "Actb"]);
            gene=:gene, coord_system="fov_1")
        @test length(pts) == 3
        @test coord_system(pts) == "fov_1"
        @test length(features(pts)) == 2
        @test "Actb" in features(pts)
        @test "Gapdh" in features(pts)
        # Actb entries share feature_id, Gapdh has a different one
        @test features(pts)[feature_ids(pts)[1]] == "Actb"
        @test feature_ids(pts)[1] == feature_ids(pts)[3]
        @test feature_ids(pts)[2] != feature_ids(pts)[1]
    end

    @testset "SpatialPoints — GeoInterface" begin
        coords = [Point2f(0, 0), Point2f(1, 0), Point2f(0, 1)]
        pts = SpatialPoints(coords)
        @test GeoInterface.isgeometry(pts)
        @test GeoInterface.geomtrait(pts) isa GeoInterface.MultiPointTrait
        @test GeoInterface.ngeom(GeoInterface.geomtrait(pts), pts) == 3
        @test GeoInterface.getgeom(GeoInterface.geomtrait(pts), pts, 1) == Point2f(0, 0)
    end

    @testset "SpatialPoints — apply (copy)" begin
        pts = SpatialPoints([Point2f(0, 0), Point2f(1, 0)]; coord_system="fov_1")
        t = translation(10.0, 20.0, "fov_1", "global")
        pts2 = apply(t, pts)
        @test pts2 isa SpatialPoints{Float32}
        @test coord_system(pts2) == "global"
        @test coords(pts2)[1] ≈ Point2f(10, 20)
        @test coords(pts2)[2] ≈ Point2f(11, 20)
        @test coords(pts)[1] == Point2f(0, 0)        # original unchanged
        @test coord_system(pts) == "fov_1"
    end

    @testset "SpatialPoints — apply! (in-place)" begin
        pts = SpatialPoints([Point2f(0, 0), Point2f(1, 0)]; coord_system="fov_1")
        t = translation(10.0, 20.0, "fov_1", "global")
        result = apply!(t, pts)
        @test result === pts                          # same object
        @test coord_system(pts) == "global"
        @test coords(pts)[1] ≈ Point2f(10, 20)
        @test coords(pts)[2] ≈ Point2f(11, 20)
    end

    @testset "SpatialShapes — construction + bbox" begin
        shp = SpatialShapes(
            [Polygon([Point2f(0,0), Point2f(1,0), Point2f(1,1), Point2f(0,1), Point2f(0,0)]),
             Polygon([Point2f(2,2), Point2f(3,2), Point2f(3,3), Point2f(2,3), Point2f(2,2)])];
            coord_system="global")
        @test shp isa SpatialShapes{<:Polygon}
        @test length(shp) == 2
        @test coord_system(shp) == "global"
        @test bbox(shp)[1, 1] ≈ 0.0   # row 1: xmin
        @test bbox(shp)[1, 2] ≈ 1.0   # row 1: xmax
        @test bbox(shp)[1, 3] ≈ 0.0   # row 1: ymin
        @test bbox(shp)[1, 4] ≈ 1.0   # row 1: ymax
        @test bbox(shp)[2, 1] ≈ 2.0   # row 2: xmin
        @test bbox(shp)[2, 2] ≈ 3.0   # row 2: xmax
    end

    @testset "SpatialShapes — GeoInterface" begin
        ring = [Point2f(0,0), Point2f(1,0), Point2f(1,1), Point2f(0,1), Point2f(0,0)]
        shp = SpatialShapes([Polygon(ring)])
        @test GeoInterface.isgeometry(shp)
        @test GeoInterface.geomtrait(shp) isa GeoInterface.GeometryCollectionTrait
        @test GeoInterface.ngeom(GeoInterface.geomtrait(shp), shp) == 1
        @test GeoInterface.getgeom(GeoInterface.geomtrait(shp), shp, 1) isa Polygon
    end

    @testset "SpatialShapes — apply (copy)" begin
        ring = [Point2f(0,0), Point2f(1,0), Point2f(1,1), Point2f(0,1), Point2f(0,0)]
        shp = SpatialShapes([Polygon(ring)]; coord_system="fov_1")
        t = translation(10.0, 20.0, "fov_1", "global")
        shp2 = apply(t, shp)
        @test coord_system(shp2) == "global"
        @test bbox(shp2)[1, 1] ≈ 10.0   # xmin shifted by 10
        @test bbox(shp2)[1, 3] ≈ 20.0   # ymin shifted by 20
        @test coord_system(shp) == "fov_1"   # original unchanged
        @test bbox(shp)[1, 1] ≈ 0.0
    end

    @testset "SpatialShapes — apply! (in-place)" begin
        ring = [Point2f(0,0), Point2f(1,0), Point2f(1,1), Point2f(0,1), Point2f(0,0)]
        shp = SpatialShapes([Polygon(ring)]; coord_system="fov_1")
        t = translation(10.0, 20.0, "fov_1", "global")
        result = apply!(t, shp)
        @test result === shp
        @test coord_system(shp) == "global"
        @test bbox(shp)[1, 1] ≈ 10.0
        @test bbox(shp)[1, 3] ≈ 20.0
    end

    @testset "Typed dataset accessors" begin
        ds = SpatialDataset()
        try
            coords = [Point2f(0, 0), Point2f(1, 1)]
            pts = SpatialPoints(coords; coord_system="global")
            ds["transcripts"] = pts
            @test points(ds, "transcripts") === pts

            ring = [Point2f(0,0), Point2f(2,0), Point2f(2,2), Point2f(0,2), Point2f(0,0)]
            shp = SpatialShapes([Polygon(ring)]; coord_system="global")
            ds["cells"] = shp
            @test shapes(ds, "cells") === shp

            @test_throws ErrorException points(ds, "cells")   # wrong type
            @test_throws ErrorException shapes(ds, "transcripts")
        finally
            close(ds)
        end
    end

end

@testset "SpatialOmics M3" begin

    using GeometryBasics
    using GeoInterface

    # ── shared fixtures ────────────────────────────────────────────────────────

    pts = SpatialPoints(
        [Point2f(x, y) for x in 0f0:1f0:4f0 for y in 0f0:1f0:4f0];
        coord_system="global")                    # 5×5 grid: 25 points

    shp = SpatialShapes(
        [Polygon([Point2f(x,y), Point2f(x+1,y), Point2f(x+1,y+1),
                  Point2f(x,y+1), Point2f(x,y)])
         for x in 0f0:2f0:4f0 for y in 0f0:2f0:4f0];
        coord_system="global")                    # 3×3 = 9 unit squares

    # ── SpatialExtent ──────────────────────────────────────────────────────────

    @testset "SpatialExtent construction" begin
        ext = SpatialExtent(0, 2, 0, 2; coord_system="global")
        @test ext.xmin == 0.0
        @test ext.xmax == 2.0
        @test coord_system(ext) == "global"
    end

    # ── SpatialROI ─────────────────────────────────────────────────────────────

    @testset "SpatialROI construction" begin
        ring = [Point2f(0,0), Point2f(2,0), Point2f(2,2), Point2f(0,2), Point2f(0,0)]
        poly = Polygon(ring)
        roi = SpatialROI(poly; coord_system="global")
        @test geometry(roi) === poly
        @test coord_system(roi) == "global"
        @test roi.extent.xmin ≈ 0.0
        @test roi.extent.xmax ≈ 2.0
    end

    # ── view on SpatialPoints with SpatialExtent ───────────────────────────────

    @testset "view(pts, SpatialExtent)" begin
        ext = SpatialExtent(0, 2, 0, 2; coord_system="global")
        v = view(pts, ext)
        @test v isa SpatialElementView{<:SpatialPoints, SpatialExtent}
        @test coord_system(v) == "global"
        # x∈{0,1,2} × y∈{0,1,2} → 9 points
        @test length(v) == 9
    end

    @testset "collect(view(pts, SpatialExtent))" begin
        ext = SpatialExtent(0, 2, 0, 2; coord_system="global")
        sub = collect(view(pts, ext))
        @test sub isa SpatialPoints
        @test length(sub) == 9
        @test coord_system(sub) == "global"
        @test all(p -> p[1] <= 2.0 && p[2] <= 2.0, coords(sub))
    end

    # ── view on SpatialShapes with SpatialExtent ───────────────────────────────

    @testset "view(shp, SpatialExtent)" begin
        ext = SpatialExtent(0, 3, 0, 3; coord_system="global")
        sub = collect(view(shp, ext))
        @test sub isa SpatialShapes
        @test coord_system(sub) == "global"
        @test length(sub) >= 1
    end

    # ── overlap=:any vs :full ──────────────────────────────────────────────────

    @testset "view(shp, SpatialExtent) — overlap modes" begin
        # ext straddles corners of several unit squares: only the (2,2) square
        # sits fully inside 1..3 × 1..3; the three corner squares touch the boundary
        ext = SpatialExtent(1, 3, 1, 3; coord_system="global")
        n_any  = length(view(shp, ext))
        n_full = length(view(shp, ext; overlap=:full))
        @test n_any > n_full
        @test n_full >= 1
        sub = collect(view(shp, ext; overlap=:full))
        bb  = bbox(sub)
        @test all(bb[:, 1] .>= 1.0)   # xmin ≥ ext.xmin
        @test all(bb[:, 2] .<= 3.0)   # xmax ≤ ext.xmax
        @test all(bb[:, 3] .>= 1.0)   # ymin ≥ ext.ymin
        @test all(bb[:, 4] .<= 3.0)   # ymax ≤ ext.ymax
    end

    @testset "view — invalid overlap raises error" begin
        ext = SpatialExtent(0, 2, 0, 2; coord_system="global")
        @test_throws ErrorException view(pts, ext; overlap=:partial)
    end

    # ── view on SpatialPoints with SpatialROI (polygon) ───────────────────────

    @testset "view(pts, SpatialROI)" begin
        tri = Polygon([Point2f(0,0), Point2f(3,0), Point2f(1.5,3), Point2f(0,0)])
        roi = SpatialROI(tri; coord_system="global")
        sub = collect(view(pts, roi))
        @test sub isa SpatialPoints
        @test length(sub) <= length(pts)
        @test all(p -> p[1] <= 3.0 && p[2] <= 3.0, coords(sub))
    end

    # ── coord system mismatch ──────────────────────────────────────────────────

    @testset "coord system mismatch" begin
        ext = SpatialExtent(0, 2, 0, 2; coord_system="other")
        @test_throws ErrorException view(pts, ext)
    end

    # ── SpatialDatasetView ─────────────────────────────────────────────────────

    @testset "SpatialDatasetView" begin
        ds = SpatialDataset()
        try
            push!(ds, CoordinateSystem("global"))
            ds["transcripts"] = pts
            ds["cells"] = shp

            ext = SpatialExtent(0, 2, 0, 2; coord_system="global")
            v_ds = view(ds, ext)
            @test v_ds isa SpatialDatasetView
            @test haskey(v_ds, "transcripts")
            @test "cells" in collect(keys(v_ds))

            @test v_ds["transcripts"] isa SpatialElementView{<:SpatialPoints}
            @test points(v_ds, "transcripts") isa SpatialElementView{<:SpatialPoints}
            @test shapes(v_ds, "cells")       isa SpatialElementView{<:SpatialShapes}
            @test length(collect(points(v_ds, "transcripts"))) == 9
        finally
            close(ds)
        end
    end

end

@testset "SpatialOmics M4" begin

    using GeometryBasics
    using Random
    Random.seed!(7)

    genes = ["Actb", "Gapdh", "Vim"]
    pts = SpatialPoints(
        (x = rand(Float32, 200) .* 500f0,
         y = rand(Float32, 200) .* 500f0,
         g = [genes[rand(1:3)] for _ in 1:200]);
        x=:x, y=:y, gene=:g, coord_system="px")

    cells = SpatialShapes(
        [let cx = rand(Float32)*450f0+25f0, cy = rand(Float32)*450f0+25f0
             Polygon([Point2f(cx-10,cy-10), Point2f(cx+10,cy-10),
                      Point2f(cx+10,cy+10), Point2f(cx-10,cy+10),
                      Point2f(cx-10,cy-10)])
         end for _ in 1:15];
        instance_id=Int32.(1:15), coord_system="px")

    @testset "SpatialPoints zarr roundtrip" begin
        path = mktempdir()
        try
            ds = SpatialDataset()
            push!(ds, CoordinateSystem("px"; units=("px","px")))
            ds["pts"] = pts
            write(ds, path, SpatialDataZarr())
            close(ds)

            ds2 = read(SpatialDataZarr(), path)
            pts2 = points(ds2, "pts")
            @test length(pts2)       == length(pts)
            @test features(pts2)     == features(pts)
            @test coords(pts2)[1]    ≈  coords(pts)[1]
            @test coord_system(pts2) == coord_system(pts)
            @test pts2.instance_id   == pts.instance_id
            close(ds2)
        finally
            rm(path; recursive=true, force=true)
        end
    end

    @testset "SpatialShapes zarr roundtrip" begin
        path = mktempdir()
        try
            ds = SpatialDataset()
            push!(ds, CoordinateSystem("px"; units=("px","px")))
            ds["cells"] = cells
            write(ds, path, SpatialDataZarr())
            close(ds)

            ds2 = read(SpatialDataZarr(), path)
            cells2 = shapes(ds2, "cells")
            @test length(cells2)       == length(cells)
            @test coord_system(cells2) == coord_system(cells)
            @test cells2.instance_id   == cells.instance_id
            r1 = GeoInterface.coordinates(geometries(cells)[1])[1]
            r2 = GeoInterface.coordinates(geometries(cells2)[1])[1]
            @test length(r1) == length(r2)
            @test all(r1[i][1] ≈ r2[i][1] && r1[i][2] ≈ r2[i][2] for i in eachindex(r1))
            close(ds2)
        finally
            rm(path; recursive=true, force=true)
        end
    end

    @testset "coord_systems preserved across roundtrip" begin
        path = mktempdir()
        try
            ds = SpatialDataset()
            push!(ds, CoordinateSystem("px"; axes=(:x,:y), units=("px","px")))
            ds["pts"] = pts
            write(ds, path, SpatialDataZarr())
            close(ds)

            ds2 = read(SpatialDataZarr(), path)
            @test "px" in coord_systems(ds2)
            close(ds2)
        finally
            rm(path; recursive=true, force=true)
        end
    end

    @testset "spill threshold triggers zarr write on attach" begin
        path = mktempdir()
        try
            ds = SpatialDataset(; path, spill_threshold=0)
            ds["pts"] = pts
            @test isfile(joinpath(path, "points", "pts", "zarr.json"))
            @test isfile(joinpath(path, "points", "pts", "coords", "zarr.json"))
            close(ds)
        finally
            rm(path; recursive=true, force=true)
        end
    end

    @testset "write produces valid zarr layout" begin
        path = mktempdir()
        try
            ds = SpatialDataset()
            ds["pts"] = pts; ds["cells"] = cells
            write(ds, path, SpatialDataZarr())
            close(ds)
            @test isfile(joinpath(path, "zarr.json"))
            @test isfile(joinpath(path, "points", "pts",   "coords",      "zarr.json"))
            @test isfile(joinpath(path, "shapes", "cells", "geom_data",   "zarr.json"))
            @test isfile(joinpath(path, "shapes", "cells", "poly_offsets","zarr.json"))
        finally
            rm(path; recursive=true, force=true)
        end
    end

end

@testset "SpatialOmics M5" begin

    using Logging

    # ── construction ─────────────────────────────────────────────────────────────

    @testset "SpatialImage 2D construction" begin
        arr = rand(Float32, 64, 64)
        img = SpatialImage(arr; coord_system="px")
        @test img.axes == (:y, :x)
        @test nchannels(img) == 1
        @test size(img) == (64, 64)
        @test isempty(img.pyramid)
        @test coord_system(img) == "px"
    end

    @testset "SpatialImage 3D construction — (c,y,x)" begin
        arr = rand(Float32, 3, 128, 128)
        img = SpatialImage(arr;
                           axes=(:c, :y, :x),
                           channel_names=["DAPI", "GFP", "RFP"],
                           coord_system="global")
        @test img.axes == (:c, :y, :x)
        @test nchannels(img) == 3
        @test channel_names(img) == ["DAPI", "GFP", "RFP"]
        @test size(img) == (3, 128, 128)
        @test coord_system(img) == "global"
    end

    @testset "SpatialImage 3D construction — (y,x,c)" begin
        arr = rand(Float32, 128, 128, 3)
        img = SpatialImage(arr;
                           axes=(:y, :x, :c),
                           channel_names=["DAPI", "GFP", "RFP"],
                           coord_system="global")
        @test img.axes == (:y, :x, :c)
        @test nchannels(img) == 3
        @test size(img) == (128, 128, 3)
    end

    @testset "SpatialImage pixel_to_cs default" begin
        arr = rand(Float32, 64, 64)
        img = SpatialImage(arr; coord_system="px")
        @test img.pixel_to_cs isa Identity
        @test img.pixel_to_cs.src == "pixel"
        @test img.pixel_to_cs.dst == "px"
    end

    @testset "SpatialImage pixel_to_cs custom" begin
        t = translation(10.0, 20.0, "pixel", "global")
        arr = rand(Float32, 64, 64)
        img = SpatialImage(arr; coord_system="global", pixel_to_cs=t)
        @test img.pixel_to_cs isa Affine
        @test img.pixel_to_cs.src == "pixel"
        @test img.pixel_to_cs.dst == "global"
    end

    # ── pyramid ───────────────────────────────────────────────────────────────────

    @testset "build_pyramid! level count" begin
        arr = rand(Float32, 3, 128, 128)
        img = SpatialImage(arr)
        build_pyramid!(img, 3)
        @test length(img.pyramid) == 3
    end

    @testset "build_pyramid! spatial dims shrink — (c,y,x)" begin
        arr = rand(Float32, 3, 128, 128)
        img = SpatialImage(arr; axes=(:c, :y, :x))
        build_pyramid!(img, 2)
        @test size(img.pyramid[1], 1) == 3         # channel dim preserved
        @test size(img.pyramid[1], 2) < 128
        @test size(img.pyramid[1], 3) < 128
        @test size(img.pyramid[2], 2) < size(img.pyramid[1], 2)
    end

    @testset "build_pyramid! spatial dims shrink — (y,x,c)" begin
        arr = rand(Float32, 128, 128, 3)
        img = SpatialImage(arr; axes=(:y, :x, :c))
        build_pyramid!(img, 2)
        @test size(img.pyramid[1], 1) < 128
        @test size(img.pyramid[1], 2) < 128
        @test size(img.pyramid[1], 3) == 3         # channel dim preserved
        @test size(img.pyramid[2], 1) < size(img.pyramid[1], 1)
    end

    @testset "build_pyramid! 2D (no channel dim)" begin
        arr = rand(Float32, 64, 64)
        img = SpatialImage(arr)
        build_pyramid!(img, 2)
        @test length(img.pyramid) == 2
        @test size(img.pyramid[1], 1) < 64
        @test size(img.pyramid[1], 2) < 64
    end

    @testset "build_pyramid! replace clears old levels" begin
        arr = rand(Float32, 3, 64, 64)
        img = SpatialImage(arr)
        build_pyramid!(img, 3)
        build_pyramid!(img, 1)
        @test length(img.pyramid) == 1
    end

    # ── zarr roundtrip ────────────────────────────────────────────────────────────

    @testset "SpatialImage zarr roundtrip — data" begin
        arr = rand(Float32, 3, 64, 64)
        img = SpatialImage(arr; coord_system="px")
        path = mktempdir()
        try
            ds = SpatialDataset()
            ds["img"] = img
            with_logger(SimpleLogger(stderr, Logging.Error)) do
                write(ds, path, SpatialDataZarr())
            end
            close(ds)

            ds2 = with_logger(SimpleLogger(stderr, Logging.Error)) do
                read(SpatialDataZarr(), path)
            end
            img2 = images(ds2, "img")
            @test size(img2.data) == (3, 64, 64)
            @test img2.data ≈ arr
            close(ds2)
        finally
            rm(path; recursive=true, force=true)
        end
    end

    @testset "SpatialImage zarr roundtrip — metadata preserved" begin
        t = translation(5.0, 10.0, "pixel", "global")
        arr = rand(Float32, 2, 32, 32)
        img = SpatialImage(arr;
                           axes=(:c, :y, :x),
                           channel_names=["ch1", "ch2"],
                           coord_system="global",
                           pixel_to_cs=t)
        path = mktempdir()
        try
            ds = SpatialDataset()
            ds["img"] = img
            with_logger(SimpleLogger(stderr, Logging.Error)) do
                write(ds, path, SpatialDataZarr())
            end
            close(ds)

            ds2 = with_logger(SimpleLogger(stderr, Logging.Error)) do
                read(SpatialDataZarr(), path)
            end
            img2 = images(ds2, "img")
            @test img2.axes == (:c, :y, :x)
            @test channel_names(img2) == ["ch1", "ch2"]
            @test coord_system(img2) == "global"
            @test img2.pixel_to_cs isa Affine
            close(ds2)
        finally
            rm(path; recursive=true, force=true)
        end
    end

    @testset "SpatialImage zarr roundtrip — pyramid preserved" begin
        arr = rand(Float32, 3, 64, 64)
        img = SpatialImage(arr; coord_system="px")
        build_pyramid!(img, 2)
        path = mktempdir()
        try
            ds = SpatialDataset()
            ds["img"] = img
            with_logger(SimpleLogger(stderr, Logging.Error)) do
                write(ds, path, SpatialDataZarr())
            end
            close(ds)

            ds2 = with_logger(SimpleLogger(stderr, Logging.Error)) do
                read(SpatialDataZarr(), path)
            end
            img2 = images(ds2, "img")
            @test length(img2.pyramid) == 2
            @test size(img2.pyramid[1]) == size(img.pyramid[1])
            close(ds2)
        finally
            rm(path; recursive=true, force=true)
        end
    end

    @testset "images accessor type error" begin
        ds = SpatialDataset()
        try
            ring = [Point2f(0,0), Point2f(1,0), Point2f(1,1), Point2f(0,1), Point2f(0,0)]
            ds["cells"] = SpatialShapes([Polygon(ring)]; coord_system="global")
            @test_throws ErrorException images(ds, "cells")
        finally
            close(ds)
        end
    end

end
