using Test
using SpatialOmicsBase
using DataFrames
using Tables
using GeometryBasics

@testset "SpatialOmicsBase.jl" begin

    @testset "SpatialDataset construction" begin
        ds = spatial_dataset()
        @test ds isa SpatialDataset{Float32}
        @test isempty(ds.images)
        @test isempty(ds.points)
        @test isempty(ds.labels)
        @test isempty(ds.shapes)
        @test isempty(ds.tables)

        ds64 = spatial_dataset(T = Float64)
        @test ds64 isa SpatialDataset{Float64}
    end

    @testset "SpatialImage" begin
        data = rand(Float32, 100, 100)
        axes = (x = 1:100, y = 1:100)
        img = SpatialImage(data, axes, Dict{String,Any}())
        @test img isa SpatialElement
        @test size(img.data) == (100, 100)
    end

    @testset "SpatialPoints" begin
        # Table constructor with label_col
        df = DataFrame(x=rand(Float32, 50), y=rand(Float32, 50),
                       gene=rand(["Spry2","Ldha","Actb"], 50), cell_id=1:50)
        pts = SpatialPoints(df; x_col=:x, y_col=:y, label_col=:gene)
        @test pts isa SpatialElement
        @test size(pts.coordinates) == (50, 2)
        @test length(pts.labels) == 50
        @test !isempty(pts.index)
        @test pts.features isa NamedTuple
        @test haskey(pts.features, :cell_id)

        # String key lookup
        sub = pts["Spry2"]
        @test sub isa SpatialPoints
        @test all(sub.labels .== "Spry2")

        # String vector lookup
        sub2 = pts[["Spry2","Ldha"]]
        @test sub2 isa SpatialPoints
        @test all(l -> l ∈ ("Spry2","Ldha"), sub2.labels)

        # Row accessor
        pt = pts[1]
        @test pt isa SpatialPoint
        @test length(pt.coordinates) == 2
        @test pt.label isa String

        # Bool mask
        mask = pts.labels .== "Actb"
        sub3 = pts[mask]
        @test all(sub3.labels .== "Actb")

        # No-label constructor
        pts_nolabel = SpatialPoints(df; x_col=:x, y_col=:y)
        @test isempty(pts_nolabel.labels)
        @test isempty(pts_nolabel.index)
    end

    @testset "setindex! routing" begin
        ds = spatial_dataset()
        df = DataFrame(x=rand(Float32, 10), y=rand(Float32, 10),
                       gene=fill("ACTB", 10))
        pts = SpatialPoints(df; x_col=:x, y_col=:y, label_col=:gene)
        ds["transcripts"] = pts
        @test haskey(ds.points, "transcripts")
        @test points(ds, "transcripts") === pts
    end

    @testset "CoordinateSystem" begin
        cs = CoordinateSystem("global")
        @test cs.name == "global"
        @test length(cs.axes) == 2

        cs3d = CoordinateSystem("3d_space", ["x","y","z"], ["µm","µm","µm"])
        @test length(cs3d.axes) == 3
    end

    @testset "IdentityTransformation" begin
        t = IdentityTransformation()
        coords = rand(Float32, 10, 2)
        out = transform_coordinates(t, coords)
        @test out == coords
        @test invert_transformation(t) isa IdentityTransformation
    end

    @testset "Chunking heuristic" begin
        sz = (1024, 1024)
        chunk = SpatialOmicsBase.default_chunk_size(sz, 4)
        @test length(chunk) == 2
        @test all(c -> c > 0, chunk)
        @test all(i -> chunk[i] <= sz[i], 1:2)
    end

    @testset "Metadata helpers" begin
        ds = spatial_dataset()
        ds.metadata["platform"] = "xenium"
        @test SpatialOmicsBase.metadata(ds)["platform"] == "xenium"
    end

    @testset "SpatialShapes" begin
        ring = Point2f[(0,0),(1,0),(1,1),(0,1),(0,0)]
        poly = GeometryBasics.Polygon(ring)

        # Direct construction — name + metadata, no type parameter on data
        s1 = SpatialShape(poly, "cell_A", (cell_id = Int32(1),))
        @test s1.geometry isa GeometryBasics.Polygon
        @test s1.name == "cell_A"
        @test s1.metadata.cell_id == Int32(1)

        shapes = SpatialShapes([s1], Dict{String,Any}())
        @test length(shapes) == 1
        @test geometry(shapes, 1) isa GeometryBasics.Polygon

        # Tables.jl interface — name is always first column
        @test Tables.istable(typeof(shapes))
        cols = Tables.columns(shapes)
        @test cols.name == ["cell_A"]
        @test cols.cell_id == [Int32(1)]

        # DataFrame conversion
        df = DataFrame(shapes)
        @test nrow(df) == 1
        @test "name" in names(df)
        @test "cell_id" in names(df)

        # Backward-compat constructor with DataFrame — name picked from :cell column
        feats = DataFrame(cell = ["B", "C"], cell_id = Int32[2, 3])
        poly2 = GeometryBasics.Polygon(Point2f[(2,0),(3,0),(3,1),(2,1),(2,0)])
        ss2 = SpatialShapes([poly, poly2], feats, Dict{String,Any}())
        @test length(ss2) == 2
        @test ss2.shapes[1].name == "B"
        @test Tables.getcolumn(ss2, :cell) == ["B", "C"]  # metadata carries original column

        # geometry/geometries accessors
        gs = geometries(ss2)
        @test length(gs) == 2
        @test gs[1] isa GeometryBasics.Polygon
    end

    @testset "add_roi! / extent(SpatialShapes)" begin
        ds  = spatial_dataset()
        add_roi!(ds, "test_roi", SpatialExtent(0.0, 10.0, 0.0, 10.0))
        @test haskey(shapes(ds), "test_roi")
        shp = shapes(ds, "test_roi")
        @test length(shp) == 1
        ext = extent(shp, 1)
        @test ext.xmin ≈ 0.0 && ext.xmax ≈ 10.0
        @test ext.ymin ≈ 0.0 && ext.ymax ≈ 10.0
    end

end
