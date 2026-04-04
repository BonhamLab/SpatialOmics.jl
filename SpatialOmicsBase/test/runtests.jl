using Test
using SpatialOmicsBase
using DataFrames

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
        img = SpatialImage{Float32}(data, axes, Dict{String,Any}())
        @test img isa SpatialElement
        @test size(img.data) == (100, 100)
    end

    @testset "SpatialPoints" begin
        coords = rand(Float32, 50, 2)
        feats = DataFrame(gene = fill("GAPDH", 50), cell_id = 1:50)
        pts = SpatialPoints{Float32}(coords, feats, Dict{String,Any}())
        @test pts isa SpatialElement
        @test size(pts.coordinates) == (50, 2)
        @test nrow(pts.features) == 50
    end

    @testset "add / get element helpers" begin
        ds = spatial_dataset()
        coords = rand(Float32, 10, 2)
        feats = DataFrame(gene = fill("ACTB", 10))
        pts = SpatialPoints{Float32}(coords, feats, Dict{String,Any}())
        add_points!(ds, "transcripts", pts)
        @test haskey(ds.points, "transcripts")
        @test get_points(ds, "transcripts") === pts
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
        set_metadata!(ds, "platform", "xenium")
        @test get_metadata(ds)["platform"] == "xenium"
    end

end
