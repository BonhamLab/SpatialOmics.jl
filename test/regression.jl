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
        t = SpatialOmics.translation(10.0, -5.0, "fov_1", "global")
        @test t isa Affine
        @test t.src == "fov_1"
        @test t.dst == "global"

        s = SpatialOmics.scaling(2.0, 2.0, "px", "µm")
        @test s isa Affine

        r = SpatialOmics.rotation(π/4, "a", "b")
        @test r isa Affine

        f = SpatialOmics.flip_y("local", "global")
        @test f isa Affine
    end

    @testset "Transformations — apply" begin
        # translation
        pts = [1.0 2.0; 3.0 4.0]   # 2×2
        t = SpatialOmics.translation(10.0, 20.0, "a", "b")
        out = apply(t, pts)
        @test out ≈ [11.0 22.0; 13.0 24.0]

        # flip_y
        f = SpatialOmics.flip_y("a", "b")
        out2 = apply(f, pts)
        @test out2 ≈ [1.0 -2.0; 3.0 -4.0]

        # identity
        id = Identity("a", "b")
        @test apply(id, pts) === pts

        # compose two translations
        t1 = SpatialOmics.translation(1.0, 0.0, "a", "b")
        t2 = SpatialOmics.translation(0.0, 1.0, "b", "c")
        tc = SpatialOmics.compose(t1, t2)
        @test apply(tc, [0.0 0.0]) ≈ [1.0 1.0]
    end

    @testset "Transformations — resolve / Dijkstra" begin
        transforms = AbstractTransformation[
            SpatialOmics.translation(100.0, 200.0, "fov_1", "global"),
            SpatialOmics.scaling(0.5, 0.5, "px", "µm"),
        ]
        t = resolve(transforms, "fov_1", "global")
        @test t isa Affine

        # no path
        @test_throws ErrorException resolve(transforms, "nowhere", "global")

        # identity (same src == dst)
        t2 = resolve(transforms, "global", "global")
        @test t2 isa Identity
    end

    @testset "Sequence apply" begin
        t1 = SpatialOmics.translation(1.0, 0.0, "a", "b")
        t2 = SpatialOmics.translation(0.0, 1.0, "b", "c")
        seq = Sequence([t1, t2], "a", "c")
        out = apply(seq, [0.0 0.0])
        @test out ≈ [1.0 1.0]
    end

    @testset "apply on SVector / Point2f" begin
        using StaticArrays
        t = SpatialOmics.translation(10.0, 20.0, "a", "b")

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
        t1 = SpatialOmics.translation(1.0, 0.0, "a", "b")
        t2 = SpatialOmics.translation(0.0, 1.0, "b", "c")
        seq = Sequence([t1, t2], "a", "c")
        svec_pts = [SVector(0.0, 0.0)]
        @test apply(seq, svec_pts)[1] ≈ SVector(1.0, 1.0)

        # flip_y on SVector
        f = SpatialOmics.flip_y("a", "b")
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
        close(ds; discard=true)
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
            close(ds; discard=true)               # should be a no-op (owned=false)
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

            t = SpatialOmics.translation(500.0, 300.0, "fov_1", "global")
            push!(ds, t)
            resolved = transform(ds, "fov_1", "global")
            @test resolved isa Affine

            pts = [0.0 0.0; 10.0 20.0]
            out = apply(resolved, pts)
            @test out ≈ [500.0 300.0; 510.0 320.0]
        finally
            close(ds; discard=true)
        end
    end

    @testset "SpatialDataset — element setindex/getindex" begin
        ds = SpatialDataset()
        try
            @test_throws MethodError setindex!(ds, (x = 1, y = 2), "test")
        finally
            close(ds; discard=true)
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
        t = SpatialOmics.translation(10.0, 20.0, "fov_1", "global")
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
        t = SpatialOmics.translation(10.0, 20.0, "fov_1", "global")
        result = apply!(t, pts)
        @test result === pts                          # same object
        @test coord_system(pts) == "global"
        @test coords(pts)[1] ≈ Point2f(10, 20)
        @test coords(pts)[2] ≈ Point2f(11, 20)
    end

    @testset "SpatialShapes — construction" begin
        shp = SpatialShapes(
            [Polygon([Point2f(0,0), Point2f(1,0), Point2f(1,1), Point2f(0,1), Point2f(0,0)]),
             Polygon([Point2f(2,2), Point2f(3,2), Point2f(3,3), Point2f(2,3), Point2f(2,2)])];
            coord_system="global")
        @test shp isa SpatialShapes{<:Polygon}
        @test length(shp) == 2
        @test coord_system(shp) == "global"
        @test length(geometries(shp)) == 2
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
        t = SpatialOmics.translation(10.0, 20.0, "fov_1", "global")
        shp2 = apply(t, shp)
        @test coord_system(shp2) == "global"
        @test GeoInterface.coordinates(geometries(shp2)[1])[1][1][1] ≈ 10.0  # x shifted
        @test GeoInterface.coordinates(geometries(shp2)[1])[1][1][2] ≈ 20.0  # y shifted
        @test coord_system(shp) == "fov_1"   # original unchanged
    end

    @testset "SpatialShapes — iteration and filter" begin
        rings = [[Point2f(i,0), Point2f(i+1,0), Point2f(i+1,1), Point2f(i,1), Point2f(i,0)]
                 for i in 0:2]
        shp = SpatialShapes(Polygon.(rings);
                            instance_id=Int32[10, 20, 30], coord_system="global")

        @test length(collect(shp)) == 3
        @test eltype(shp) <: SpatialShape

        row = shp[2]
        @test row isa SpatialShape
        @test row.instance_id == Int32(20)
        @test row.coord_system == "global"

        kept = filter(s -> s.instance_id in [10, 30], shp)
        @test length(kept) == 2
        @test kept.instance_id == Int32[10, 30]
        @test coord_system(kept) == "global"
    end

    @testset "SpatialShapes — apply! (in-place)" begin
        ring = [Point2f(0,0), Point2f(1,0), Point2f(1,1), Point2f(0,1), Point2f(0,0)]
        shp = SpatialShapes([Polygon(ring)]; coord_system="fov_1")
        t = SpatialOmics.translation(10.0, 20.0, "fov_1", "global")
        result = apply!(t, shp)
        @test result === shp
        @test coord_system(shp) == "global"
        @test GeoInterface.coordinates(geometries(shp)[1])[1][1][1] ≈ 10.0
        @test GeoInterface.coordinates(geometries(shp)[1])[1][1][2] ≈ 20.0
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
            close(ds; discard=true)
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

    @testset "SpatialExtent from SpatialShapes" begin
        rings = [[Point2f(0,0), Point2f(1,0), Point2f(1,1), Point2f(0,1), Point2f(0,0)],
                 [Point2f(2,2), Point2f(4,2), Point2f(4,5), Point2f(2,5), Point2f(2,2)]]
        shp_ext = SpatialShapes(Polygon.(rings); instance_id=Int32[1,2], coord_system="global")
        ext = SpatialExtent(shp_ext)
        @test ext.xmin == 0.0 && ext.xmax == 4.0
        @test ext.ymin == 0.0 && ext.ymax == 5.0
        @test coord_system(ext) == "global"

        # filter then extent — the idiomatic pipeline
        ext2 = SpatialExtent(filter(s -> s.instance_id == Int32(1), shp_ext))
        @test ext2.xmax == 1.0
    end

    @testset "SpatialExtent union" begin
        a = SpatialExtent(0, 2, 0, 2; coord_system="g")
        b = SpatialExtent(1, 4, 1, 3; coord_system="g")
        u = a ∪ b
        @test u.xmin == 0.0 && u.xmax == 4.0
        @test u.ymin == 0.0 && u.ymax == 3.0
        @test coord_system(u) == "g"
        @test_throws ErrorException SpatialExtent(0,1,0,1;coord_system="a") ∪
                                    SpatialExtent(0,1,0,1;coord_system="b")
    end

    @testset "SpatialExtent intersect" begin
        a = SpatialExtent(0, 3, 0, 3; coord_system="g")
        b = SpatialExtent(1, 4, 1, 4; coord_system="g")
        i = a ∩ b
        @test i isa SpatialExtent
        @test i.xmin == 1.0 && i.xmax == 3.0
        @test i.ymin == 1.0 && i.ymax == 3.0
        # non-overlapping
        c = SpatialExtent(5, 6, 5, 6; coord_system="g")
        @test isnothing(a ∩ c)
        # touching at edge — not an overlap
        d = SpatialExtent(3, 5, 0, 3; coord_system="g")
        @test isnothing(a ∩ d)
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
        # all coords of fully-inside shapes lie within [1,3]×[1,3]
        @test all(geometries(sub)) do g
            all(Iterators.flatten(GeoInterface.coordinates(g))) do pt
                1.0 <= pt[1] <= 3.0 && 1.0 <= pt[2] <= 3.0
            end
        end
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
            close(ds; discard=true)
        end
    end

    # ── SpatialShapes(::SpatialExtent) and SpatialShapes(::SpatialROI) ──────────

    @testset "SpatialShapes(SpatialExtent)" begin
        ext = SpatialExtent(1.0, 3.0, 2.0, 5.0; coord_system="global")
        s   = SpatialShapes(ext)
        @test s isa SpatialShapes
        @test length(s) == 1
        @test s.instance_id == Int32[1]
        @test coord_system(s) == "global"
        ring = GeoInterface.coordinates(s.geometries[1])[1]
        @test length(ring) == 5
        @test ring[1] ≈ ring[end]
        xs = [p[1] for p in ring]; ys = [p[2] for p in ring]
        @test minimum(xs) ≈ 1.0 && maximum(xs) ≈ 3.0
        @test minimum(ys) ≈ 2.0 && maximum(ys) ≈ 5.0
    end

    @testset "SpatialShapes(SpatialROI)" begin
        ring = [Point2f(0,0), Point2f(2,0), Point2f(1,2), Point2f(0,0)]
        roi  = SpatialROI(Polygon(ring); coord_system="global")
        s    = SpatialShapes(roi)
        @test length(s) == 1
        @test coord_system(s) == "global"
        r2 = GeoInterface.coordinates(s.geometries[1])[1]
        @test r2[1] ≈ r2[end]
    end

    @testset "SpatialShapes(SpatialExtent) zarr roundtrip" begin
        ext  = SpatialExtent(0.0, 10.0, 0.0, 10.0; coord_system="global")
        path = mktempdir()
        try
            ds = SpatialDataset()
            push!(ds, CoordinateSystem("global"))
            ds["roi"] = SpatialShapes(ext)
            write!(ds, path, SpatialDataZarr())
            close(ds; discard=true)
            ds2 = read(SpatialDataZarr(), path)
            s2  = shapes(ds2, "roi")
            @test length(s2) == 1
            @test coord_system(s2) == "global"
            close(ds2; discard=true)
        finally
            rm(path; recursive=true, force=true)
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
        mktempdir() do path
            ds = SpatialDataset()
            push!(ds, CoordinateSystem("px"; units=("px","px")))
            ds["pts"] = pts
            write(ds, path, SpatialDataZarr())
            close(ds; discard=true)

            ds2 = read(SpatialDataZarr(), path)
            pts2 = points(ds2, "pts")
            @test length(pts2)       == length(pts)
            @test features(pts2)     == features(pts)
            @test coords(pts2)[1]    ≈  coords(pts)[1]
            @test coord_system(pts2) == coord_system(pts)
            @test pts2.instance_id   == pts.instance_id
            close(ds2; discard=true)
        end
    end

    @testset "SpatialShapes zarr roundtrip" begin
        mktempdir() do path
            ds = SpatialDataset()
            push!(ds, CoordinateSystem("px"; units=("px","px")))
            ds["cells"] = cells
            write(ds, path, SpatialDataZarr())
            close(ds; discard=true)

            ds2 = read(SpatialDataZarr(), path)
            cells2 = shapes(ds2, "cells")
            @test length(cells2)       == length(cells)
            @test coord_system(cells2) == coord_system(cells)
            @test cells2.instance_id   == cells.instance_id
            r1 = GeoInterface.coordinates(geometries(cells)[1])[1]
            r2 = GeoInterface.coordinates(geometries(cells2)[1])[1]
            @test length(r1) == length(r2)
            @test all(r1[i][1] ≈ r2[i][1] && r1[i][2] ≈ r2[i][2] for i in eachindex(r1))
            close(ds2; discard=true)
        end
    end

    @testset "coord_systems preserved across roundtrip" begin
        mktempdir() do path
            ds = SpatialDataset(; path)
            push!(ds, CoordinateSystem("px"; axes=(:x,:y), units=("px","px")))
            @test isdirty(ds)
            save!(ds)
            ds2 = read(SpatialDataZarr(), path)
            @test "px" in coord_systems(ds2)
            close(ds; discard=true); close(ds2; discard=true)
        end
    end

    @testset "transforms preserved across roundtrip" begin
        mktempdir() do path
            ds = SpatialDataset(; path)
            push!(ds, CoordinateSystem("fov"; axes=(:x,:y), units=("µm","µm")))
            push!(ds, CoordinateSystem("global"; axes=(:x,:y), units=("µm","µm")))
            push!(ds, SpatialOmics.translation(100.0, 200.0, "fov", "global"))
            save!(ds)
            ds2 = read(SpatialDataZarr(), path)
            @test length(ds2.transforms) == 1
            t = ds2.transforms[1]
            @test t isa Affine
            @test t.src == "fov" && t.dst == "global"
            @test transform(ds2, "fov", "global") isa AbstractTransformation
            close(ds; discard=true); close(ds2; discard=true)
        end
    end

    @testset "setindex! stages until save!" begin
        mktempdir() do path
            ds = SpatialDataset(; path)
            ds["pts"] = copy(pts)
            @test !isfile(joinpath(path, "points", "pts", "zarr.json"))
            @test isdirty(ds)
            save!(ds)
            @test isfile(joinpath(path, "points", "pts", "zarr.json"))
            @test isfile(joinpath(path, "points", "pts", "coords", "zarr.json"))
            close(ds; discard=true)
        end
    end

    @testset "metadata NamedTuple-of-vectors roundtrip" begin
        mktempdir() do path
            ds = SpatialDataset(; path)
            ds.metadata["ann"] = (fov=Int32[1, 1, 2], z=Float32[0.5, 1.0, 0.5],
                                   comp=["Cytoplasm", "Nucleus", "Cytoplasm"])
            @test !isdir(joinpath(path, "metadata", "ann"))
            save!(ds)
            @test isdir(joinpath(path, "metadata", "ann"))

            ds2 = read(SpatialDataZarr(), path)
            ann = ds2.metadata["ann"]
            @test ann.fov  == Int32[1, 1, 2]
            @test ann.z    ≈  Float32[0.5, 1.0, 0.5]
            @test ann.comp == ["Cytoplasm", "Nucleus", "Cytoplasm"]
            close(ds; discard=true)
            close(ds2; discard=true)
        end
    end

    @testset "write produces valid zarr layout" begin
        mktempdir() do path
            ds = SpatialDataset()
            ds["pts"] = copy(pts); ds["cells"] = copy(cells)
            write(ds, path, SpatialDataZarr())
            close(ds; discard=true)
            @test isfile(joinpath(path, "zarr.json"))
            @test isfile(joinpath(path, "points", "pts",   "coords",      "zarr.json"))
            @test isfile(joinpath(path, "shapes", "cells", "geom_data",   "zarr.json"))
            @test isfile(joinpath(path, "shapes", "cells", "poly_offsets","zarr.json"))
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
        t = SpatialOmics.translation(10.0, 20.0, "pixel", "global")
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
            close(ds; discard=true)

            ds2 = with_logger(SimpleLogger(stderr, Logging.Error)) do
                read(SpatialDataZarr(), path)
            end
            img2 = images(ds2, "img")
            @test size(img2.data) == (3, 64, 64)
            @test img2.data ≈ arr
            close(ds2; discard=true)
        finally
            rm(path; recursive=true, force=true)
        end
    end

    @testset "SpatialImage zarr roundtrip — metadata preserved" begin
        t = SpatialOmics.translation(5.0, 10.0, "pixel", "global")
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
            close(ds; discard=true)

            ds2 = with_logger(SimpleLogger(stderr, Logging.Error)) do
                read(SpatialDataZarr(), path)
            end
            img2 = images(ds2, "img")
            @test img2.axes == (:c, :y, :x)
            @test channel_names(img2) == ["ch1", "ch2"]
            @test coord_system(img2) == "global"
            @test img2.pixel_to_cs isa Affine
            close(ds2; discard=true)
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
            close(ds; discard=true)

            ds2 = with_logger(SimpleLogger(stderr, Logging.Error)) do
                read(SpatialDataZarr(), path)
            end
            img2 = images(ds2, "img")
            @test length(img2.pyramid) == 2
            @test size(img2.pyramid[1]) == size(img.pyramid[1])
            close(ds2; discard=true)
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
            close(ds; discard=true)
        end
    end

end

@testset "SpatialOmics M6" begin

    using Logging

    # ── shared fixtures ────────────────────────────────────────────────────────
    Random.seed!(42)
    n_cells = 50
    n_genes = 4
    genes   = ["Actb", "Gapdh", "Col1a1", "Vim"]

    cells = SpatialShapes(
        [let cx = Float32(rand()*800+100), cy = Float32(rand()*800+100)
             Polygon([Point2f(cx-30,cy-30), Point2f(cx+30,cy-30),
                      Point2f(cx+30,cy+30), Point2f(cx-30,cy+30),
                      Point2f(cx-30,cy-30)])
         end for _ in 1:n_cells];
        instance_id=Int32.(1:n_cells), coord_system="global_px")

    X   = rand(Float32, n_cells, n_genes)
    rel = SpatialRelation(Expression(), "cells", Int32.(1:n_cells), X;
                          obs=(instance_id=Int32.(1:n_cells),), var=(name=genes,))

    # ── SpatialRelation construction ───────────────────────────────────────────

    @testset "SpatialRelation construction" begin
        @test nobs(rel) == n_cells
        @test nvar(rel) == n_genes
        @test var_names(rel) == genes
        @test rel.src == "cells"
        @test rel.kind isa Expression
    end

    @testset "SpatialRelation show" begin
        s = sprint(show, rel)
        @test contains(s, "SpatialRelation")
        @test contains(s, string(n_cells))
        @test contains(s, "cells")
    end

    # ── expression weight lookup ───────────────────────────────────────────────

    @testset "expression weight lookup" begin
        actb_col = findfirst(==("Actb"), genes)
        @test rel[:, "Actb"] ≈ X[:, actb_col]
        @test rel[1, "Actb"] isa Float32
        @test rel[[1, 2], "Actb"] isa Vector{Float32}
    end

    @testset "var_names lookup" begin
        @test var_names(rel) == genes
        @test length(var_names(rel)) == n_genes
    end

    # ── passthrough accessors on SpatialElementView ────────────────────────────

    @testset "geometries on SpatialElementView" begin
        ext  = SpatialExtent(0, 500, 0, 500; coord_system="global_px")
        v    = view(cells, ext)
        geoms = geometries(v)
        @test length(geoms) == length(v)
        @test geoms isa Vector
    end

    @testset "instance_id on SpatialElementView" begin
        ext = SpatialExtent(0, 500, 0, 500; coord_system="global_px")
        v   = view(cells, ext)
        ids = instance_id(v)
        @test length(ids) == length(v)
        @test ids isa Vector{Int32}
    end

    # ── SpatialLabels ──────────────────────────────────────────────────────────

    @testset "SpatialLabels construction" begin
        data = zeros(Int32, 64, 64)
        data[10:30, 10:30] .= 1
        data[40:60, 40:60] .= 2
        lbl  = SpatialLabels(data;
                              instance_map=Dict{Int32,Int32}(1=>1, 2=>2),
                              coord_system="global_px")
        @test lbl.axes == (:y, :x)
        @test coord_system(lbl) == "global_px"
        @test length(instance_ids(lbl)) == 2
        @test size(lbl) == (64, 64)
    end

    # ── Zarr round-trip ────────────────────────────────────────────────────────

    @testset "SpatialRelation zarr roundtrip" begin
        path = mktempdir()
        try
            ds = SpatialDataset()
            ds["cells"] = cells
            ds["expr"]  = rel
            with_logger(SimpleLogger(stderr, Logging.Error)) do
                write(ds, path, SpatialDataZarr())
            end
            close(ds; discard=true)

            ds2 = with_logger(SimpleLogger(stderr, Logging.Error)) do
                read(SpatialDataZarr(), path)
            end
            rel2 = relations(ds2, "expr")
            @test nobs(rel2) == n_cells
            @test nvar(rel2) == n_genes
            @test var_names(rel2) == genes
            @test rel2.src == "cells"
            @test rel2.weights ≈ X   atol=1e-5
            close(ds2; discard=true)
        finally
            rm(path; recursive=true, force=true)
        end
    end

    @testset "SpatialLabels zarr roundtrip" begin
        path = mktempdir()
        try
            data = rand(Int32.(0:5), 32, 32)
            lbl  = SpatialLabels(data;
                                  instance_map=Dict{Int32,Int32}(i=>i for i in 1:5),
                                  coord_system="global_px")
            ds = SpatialDataset(); ds["seg"] = lbl
            with_logger(SimpleLogger(stderr, Logging.Error)) do
                write(ds, path, SpatialDataZarr())
            end
            close(ds; discard=true)

            ds2 = with_logger(SimpleLogger(stderr, Logging.Error)) do
                read(SpatialDataZarr(), path)
            end
            lbl2 = labels(ds2, "seg")
            @test size(lbl2.data) == (32, 32)
            @test lbl2.data == data
            @test coord_system(lbl2) == "global_px"
            close(ds2; discard=true)
        finally
            rm(path; recursive=true, force=true)
        end
    end

end  # M6

@testset "SpatialOmics M7" begin

    using Random, GeometryBasics, Logging

    @testset "Round-trip — all element kinds" begin
        path = mktempdir(; prefix="so_m7_roundtrip_")
        try
            Random.seed!(1)
            genes = ["Actb", "Gapdh", "Col1a1"]
            n = 20

            pts = SpatialPoints(
                [Point2f(rand()*100, rand()*100) for _ in 1:n];
                feature_id       = Int32.(rand(1:3, n)),
                feature_codebook = genes,
                instance_id      = zeros(Int32, n),
                coord_system     = "global")

            polys = [let cx=rand()*80+10f0, cy=rand()*80+10f0
                         Polygon([Point2f(cx-5,cy-5), Point2f(cx+5,cy-5),
                                  Point2f(cx+5,cy+5), Point2f(cx-5,cy+5),
                                  Point2f(cx-5,cy-5)])
                     end for _ in 1:10]
            shp = SpatialShapes(polys; instance_id=Int32.(1:10), coord_system="global")

            X   = rand(Float32, 10, 3)
            tbl = SpatialRelation(Expression(), "cells", Int32.(1:10), X;
                                  obs=(instance_id=Int32.(1:10),), var=(name=genes,))

            img = SpatialImage(rand(UInt16, 8, 8, 2);
                axes=(:y,:x,:c), channel_names=["DAPI","GFP"], coord_system="global")

            data_lbl = Int32.(rand(0:5, 8, 8))
            lbl = SpatialLabels(data_lbl;
                instance_map=Dict{Int32,Int32}(i=>i for i in 1:5),
                coord_system="global")

            ds = SpatialDataset()
            push!(ds, CoordinateSystem("global"; units=("µm","µm")))
            ds["transcripts"] = pts
            ds["cells"]       = shp
            ds["expression"]  = tbl
            ds["dapi"]        = img
            ds["seg"]         = lbl

            with_logger(SimpleLogger(stderr, Logging.Error)) do
                write(ds, path, SpatialDataZarr())
            end
            close(ds; discard=true)

            ds2 = with_logger(SimpleLogger(stderr, Logging.Error)) do
                read(SpatialDataZarr(), path)
            end

            @test haskey(ds2.elements, "transcripts")
            @test haskey(ds2.elements, "cells")
            @test haskey(ds2.relations, "expression")
            @test haskey(ds2.elements, "dapi")
            @test haskey(ds2.elements, "seg")

            pts2 = points(ds2, "transcripts")
            @test length(pts2) == n
            @test features(pts2) == genes

            shp2 = shapes(ds2, "cells")
            @test length(shp2) == 10

            tbl2 = relations(ds2, "expression")
            @test nobs(tbl2) == 10
            @test nvar(tbl2) == 3
            @test var_names(tbl2) == genes
            @test tbl2.weights ≈ X

            img2 = images(ds2, "dapi")
            @test nchannels(img2) == 2
            @test channel_names(img2) == ["DAPI","GFP"]

            lbl2 = labels(ds2, "seg")
            @test size(lbl2.data) == (8, 8)
            @test lbl2.data == data_lbl

            close(ds2; discard=true)
        finally
            rm(path; recursive=true, force=true)
        end
    end

    xenium_path = "/home/kevin/Repos/stx_dev/test_data/experiments/xenium_ex.zarr"
    if isdir(xenium_path)
        @testset "Python SpatialData read — Xenium smoke test" begin
            ds = with_logger(SimpleLogger(stderr, Logging.Error)) do
                read(SpatialDataZarr(), xenium_path)
            end

            @test haskey(ds.elements, "morphology_focus")
            @test haskey(ds.elements, "cell_labels")
            @test haskey(ds.elements, "cell_boundaries")
            @test haskey(ds.elements, "transcripts")
            @test haskey(ds.relations, "table")

            img = images(ds, "morphology_focus")
            @test img isa SpatialImage
            @test nchannels(img) == 4
            @test length(img.pyramid) >= 1

            shp = shapes(ds, "cell_boundaries")
            @test length(shp) > 0

            pts = points(ds, "transcripts")
            @test length(pts) > 0
            @test length(features(pts)) > 0

            tbl = relations(ds, "table")
            @test nvar(tbl) == 377
            @test nobs(tbl) > 0
            @test length(var_names(tbl)) == 377
            @test all(!isempty, var_names(tbl))
            @test length(obs_names(tbl)) == nobs(tbl)
        end
    end

end

@testset "SpatialOmics M8" begin

    cosmx_path = "/home/kevin/Repos/stx_dev/test_data/experiments/cosmx_ex_raw/flatFiles/mw_mus_p1_11"
    if isdir(cosmx_path)
        @testset "CosMx reader — smoke test" begin
            ds = read(CosMx(), cosmx_path)

            @test haskey(ds.elements, "transcripts")
            @test haskey(ds.elements, "cells")
            @test haskey(ds.elements, "fovs")

            pts = points(ds, "transcripts")
            @test length(pts) > 0
            @test length(features(pts)) > 0
            @test coord_system(pts) == "global_px"

            shp = shapes(ds, "cells")
            @test length(shp) > 0
            @test coord_system(shp) == "global_px"

            # per-FOV coord systems and transforms registered
            @test haskey(ds.coord_systems, "global_px")
            @test any(cs -> startswith(cs, "fov_"), keys(ds.coord_systems))
            @test any(t -> startswith(t.src, "fov_") && t.dst == "global_px",
                      ds.transforms)

            fovshp = shapes(ds, "fovs")
            n_fovs = count(cs -> startswith(cs, "fov_"), keys(ds.coord_systems))
            @test length(fovshp) == n_fovs
            @test coord_system(fovshp) == "global_px"
            @test length(sources(ds)) == n_fovs
            @test origin_ids(pts) !== nothing
            @test origin_ids(shp) !== nothing
            @test length(features(pts, :z)) == length(pts)
            @test length(features(pts, :CellComp)) == length(pts)
            @test !haskey(ds.metadata, "transcripts_annotations")
            first_source = first(sources(ds))
            @test length(points(view(ds, first_source), "transcripts")) > 0
        end
    end

end

@testset "SpatialOmics M11" begin

    using GeometryBasics, Random

    # ── fixtures: 3 square cells, 6 known transcripts ─────────────────────────
    cells = SpatialShapes(
        [Polygon([Point2f(0,0),  Point2f(10,0),  Point2f(10,10),  Point2f(0,10),  Point2f(0,0)]),
         Polygon([Point2f(20,0), Point2f(30,0),  Point2f(30,10),  Point2f(20,10), Point2f(20,0)]),
         Polygon([Point2f(40,0), Point2f(50,0),  Point2f(50,10),  Point2f(40,10), Point2f(40,0)])];
        instance_id=Int32.([1, 2, 3]))

    pts = SpatialPoints(
        [Point2f(5,5),  Point2f(5,5),
         Point2f(25,5), Point2f(25,5),
         Point2f(45,5), Point2f(45,5)];
        feature_id=Int32.([1,2,1,2,1,2]),
        feature_codebook=["GeneA","GeneB"])

    # ── analyze(Expression()) ─────────────────────────────────────────────────

    @testset "analyze Expression" begin
        rel = analyze(Expression(), pts, cells)
        @test rel.kind isa Expression
        @test nobs(rel) == 3
        @test nvar(rel) == 2
        @test size(rel.weights) == (3, 2)
        @test all(rel.weights .== 1f0)
        @test var_names(rel) == ["GeneA","GeneB"]
        @test length(rel.src_ids) == 3
    end

    # ── analyze(Membership()) ─────────────────────────────────────────────────

    @testset "analyze Membership" begin
        rel = analyze(Membership(), pts, cells)
        @test rel.kind isa Membership
        @test nobs(rel) == 6
        @test rel.dst_ids == Int32[1,1,2,2,3,3]
        @test rel.weights === nothing
    end

    # ── default dispatch ──────────────────────────────────────────────────────

    @testset "default dispatch pts+shapes → Expression" begin
        rel = analyze(pts, cells)
        @test rel.kind isa Expression
    end

    @testset "default dispatch shapes+shapes → Membership" begin
        rel = analyze(cells, cells)
        @test rel.kind isa Membership
    end

    # ── annotate — pure, shared weights ──────────────────────────────────────

    @testset "annotate" begin
        rel    = analyze(Expression(), pts, cells)
        labels = ["T","B","M"]
        rel2   = annotate(rel, labels; key=:cell_type)
        @test hasproperty(rel2.obs, :cell_type)
        @test rel2.obs.cell_type == labels
        @test rel2.weights === rel.weights        # no copy
    end

    # ── getindex — Expression relation ───────────────────────────────────────

    @testset "getindex — scalar and slices by instance_id" begin
        rel = analyze(Expression(), pts, cells)
        # scalar: one cell, one gene
        @test rel[1, "GeneA"] isa Float32
        @test rel[1, "GeneA"] == 1f0
        @test rel[2, "GeneB"] == 1f0
        # all cells, one gene → Vector
        v = rel[:, "GeneA"]
        @test v isa Vector{Float32}
        @test length(v) == 3
        @test all(v .== 1f0)
        # one cell, all genes → Vector
        r = rel[1, :]
        @test r isa Vector{Float32}
        @test length(r) == 2
        # multi-row by instance_ids, one gene → Vector
        v2 = rel[[1, 3], "GeneA"]
        @test v2 isa Vector{Float32}
        @test length(v2) == 2
        # multi-row, multi-gene → Matrix
        M = rel[[1, 3], ["GeneA", "GeneB"]]
        @test M isa Matrix{Float32}
        @test size(M) == (2, 2)
        # all cells, multi-gene → Matrix
        @test size(rel[:, ["GeneA", "GeneB"]]) == (3, 2)
        # multi-row, all genes → Matrix
        @test size(rel[[1, 2], :]) == (2, 2)
    end

    @testset "getindex — string row via obs.name" begin
        rel      = analyze(Expression(), pts, cells)
        rel_named = annotate(rel, ["cell_T", "cell_B", "cell_M"]; key=:name)
        @test rel_named["cell_T", "GeneA"] isa Float32
        @test rel_named["cell_T", "GeneA"] == 1f0
        v = rel_named[["cell_T", "cell_M"], :]
        @test size(v) == (2, 2)
        @test all(v .== 1f0)
    end

    @testset "getindex — error cases" begin
        rel = analyze(Expression(), pts, cells)
        @test_throws ErrorException rel[99, "GeneA"]         # instance_id not found
        @test_throws ErrorException rel[1, "NoGene"]         # gene name not found
        @test_throws ErrorException rel["cell_T", "GeneA"]   # no obs.name column
    end

    # ── obs_names ─────────────────────────────────────────────────────────────

    @testset "obs_names" begin
        rel = analyze(Expression(), pts, cells)
        # no obs.name → fallback to string.(src_ids)
        @test obs_names(rel) == string.(rel.src_ids)
        @test length(obs_names(rel)) == 3
        # with obs.name via annotate
        rel_named = annotate(rel, ["T_cell", "B_cell", "Mac"]; key=:name)
        @test obs_names(rel_named) == ["T_cell", "B_cell", "Mac"]
    end

    # ── multi-level analyze ───────────────────────────────────────────────────

    @testset "multi-level analyze: cells → ROIs by cell_type" begin
        # roi1 covers cells 1+2 (centroids at (5,5) and (25,5))
        # roi2 covers cell 3  (centroid at (45,5))
        rois = SpatialShapes(
            [Polygon([Point2f(-1,-1), Point2f(35,-1), Point2f(35,11),
                      Point2f(-1,11), Point2f(-1,-1)]),
             Polygon([Point2f(35,-1), Point2f(55,-1), Point2f(55,11),
                      Point2f(35,11), Point2f(35,-1)])];
            instance_id=Int32.([10, 20]))
        cell_obs = (cell_type = ["TypeA", "TypeA", "TypeB"],)
        rel2 = analyze(cells, rois, cell_obs; by=:cell_type)
        @test rel2.kind isa Expression
        @test nobs(rel2) == 2
        @test nvar(rel2) == 2
        @test sort(var_names(rel2)) == ["TypeA", "TypeB"]
        typeA_col = findfirst(==("TypeA"), var_names(rel2))
        typeB_col = findfirst(==("TypeB"), var_names(rel2))
        roi1_row  = findfirst(==(Int32(10)), rel2.src_ids)
        roi2_row  = findfirst(==(Int32(20)), rel2.src_ids)
        @test rel2.weights[roi1_row, typeA_col] == 2f0
        @test rel2.weights[roi1_row, typeB_col] == 0f0
        @test rel2.weights[roi2_row, typeA_col] == 0f0
        @test rel2.weights[roi2_row, typeB_col] == 1f0
    end

    # ── distances ─────────────────────────────────────────────────────────────

    @testset "distances shapes→shapes" begin
        d = distances(cells, cells)
        @test length(d) == length(cells)
        @test all(d .>= 0f0)
    end

    # ── PointDensity and ShapeColorView struct construction ───────────────────

    @testset "PointDensity construction" begin
        pd = density(pts; resolution=64, feature="GeneA")
        @test pd isa PointDensity
        @test pd.resolution == 64
        @test pd.feature == "GeneA"
    end

    @testset "ShapeColorView construction" begin
        rel  = analyze(Expression(), pts, cells)
        rel2 = annotate(rel, ["T","B","M"]; key=:cell_type)
        scv  = ShapeColorView(cells, rel2, :cell_type, :tab10)
        @test scv isa ShapeColorView
        @test scv.color_by == :cell_type
    end

    # ── SpatialRelation zarr round-trip ───────────────────────────────────────

    @testset "analyze + zarr round-trip" begin
        path = mktempdir()
        try
            ds   = SpatialDataset()
            ds["cells"] = cells   # attach before analyze so _element_name resolves
            rel  = analyze(Expression(), pts, cells)
            rel  = annotate(rel, ["T","B","M"]; key=:cell_type)
            ds["expr"]  = rel
            write(ds, path, SpatialDataZarr())
            close(ds; discard=true)
            ds2  = read(SpatialDataZarr(), path)
            rel2 = relations(ds2, "expr")
            @test nobs(rel2) == 3
            @test nvar(rel2) == 2
            @test rel2.weights ≈ rel.weights  atol=1e-5
            @test rel2.src == "cells"
            close(ds2; discard=true)
        finally
            rm(path; recursive=true, force=true)
        end
    end

end  # M11

# ── M12 — Real data fixtures ──────────────────────────────────────────────────
# These tests use committed zarr fixtures in test/data/ and run unconditionally
# in CI. Generate the fixtures locally with: julia --project=. test/make_fixtures.jl

@testset "SpatialOmics M12 — Real data fixtures" begin

    xenium_path = joinpath(@__DIR__, "data", "xenium_small.zarr")
    @testset "Xenium fixture" begin
        if !isdir(xenium_path)
            @warn "Xenium fixture not found at $xenium_path — run test/make_fixtures.jl to generate it"
        else
            ds  = read(SpatialDataZarr(), xenium_path)
            tx  = points(ds, "transcripts")
            shp = shapes(ds, "cell_boundaries")
            img = images(ds, "morphology_focus")
            lbl = labels(ds, "cell_labels")

            @test length(coords(tx)) > 500
            @test length(top_features(tx, 5)) == 5
            @test length(tx.instance_id) == length(coords(tx))  # structure check; values may be 0 for Python-source fixtures

            @test length(geometries(shp)) > 10

            @test nchannels(img) == 4
            @test length(channel_names(img)) == 4

            @test size(data(lbl), 1) > 0
            @test size(data(lbl), 2) > 0

            @test !isempty(coord_systems(ds))

            close(ds; discard=true)
        end
    end

    visium_path = joinpath(@__DIR__, "data", "visium_small.zarr")
    @testset "Visium fixture" begin
        if !isdir(visium_path)
            @warn "Visium fixture not found at $visium_path — run test/make_fixtures.jl to generate it"
        else
            ds  = read(SpatialDataZarr(), visium_path)
            shp = shapes(ds, "Visium_HD_Mouse_Small_Intestine_square_016um")
            img = images(ds, "Visium_HD_Mouse_Small_Intestine_lowres_image")

            @test length(geometries(shp)) > 50

            @test ndims(data(img)) >= 2
            @test size(data(img), 1) > 0

            close(ds; discard=true)
        end
    end

end  # M12
