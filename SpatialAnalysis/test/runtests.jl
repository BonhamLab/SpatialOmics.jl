using Test
using SpatialAnalysis
using SpatialOmicsBase
using DataFrames
using GeometryBasics
using Tables

@testset "SpatialAnalysis.jl" begin

    @testset "Sub-module availability" begin
        @test isdefined(SpatialAnalysis, :QualityControl)
        @test isdefined(SpatialAnalysis, :Preprocessing)
        @test isdefined(SpatialAnalysis, :CellInteraction)
        # SpatialStatistics is not a submodule — stats functions live in SpatialAnalysis directly
        @test !isdefined(SpatialAnalysis, :SpatialStatistics)
    end

    @testset "Exported symbol availability" begin
        # QC
        @test isdefined(SpatialAnalysis, :filter_spots)
        @test isdefined(SpatialAnalysis, :detect_artifacts)
        @test isdefined(SpatialAnalysis, :spatially_aware_qc)
        # Preprocessing
        @test isdefined(SpatialAnalysis, :normalize_counts)
        @test isdefined(SpatialAnalysis, :correct_batch_effects)
        @test isdefined(SpatialAnalysis, :impute_missing)
        @test isdefined(SpatialAnalysis, :crop)
        # Statistics
        @test isdefined(SpatialAnalysis, :spatial_autocorrelation)
        @test isdefined(SpatialAnalysis, :hotspot_detection)
        @test isdefined(SpatialAnalysis, :spatial_clustering)
        @test isdefined(SpatialAnalysis, :distances)
        @test isdefined(SpatialAnalysis, :classify_cells)
        # Interaction
        @test isdefined(SpatialAnalysis, :ligand_receptor_analysis)
        @test isdefined(SpatialAnalysis, :neighborhood_analysis)
        @test isdefined(SpatialAnalysis, :communication_score)
        # Pipeline
        @test isdefined(SpatialAnalysis, :run_qc_pipeline)
        @test isdefined(SpatialAnalysis, :default_qc_pipeline)
    end

    @testset "Stub functions raise NotImplemented" begin
        ds = spatial_dataset()
        @test_throws ErrorException filter_spots(ds)
        @test_throws ErrorException detect_artifacts(ds)
        @test_throws ErrorException spatially_aware_qc(ds)
        @test_throws ErrorException normalize_counts(ds)
        @test_throws ErrorException spatial_autocorrelation(ds, "GAPDH")
        @test_throws ErrorException spatial_clustering(ds)
    end

    @testset "default_qc_pipeline returns a vector of functions" begin
        steps = default_qc_pipeline()
        @test steps isa Vector{Function}
        @test length(steps) >= 1
    end

    # -------------------------------------------------------------------------
    # Helpers
    # -------------------------------------------------------------------------
    function _make_pts()
        coords = Float32[1 1; 2 2; 5 5; 6 6]
        feats  = DataFrame(gene = ["A","B","C","D"], cell_id = 1:4)
        SpatialPoints(coords, feats, Dict{String,Any}())
    end
    function _make_shp()
        geoms = [Circle(Point2f(1,1), 0.5f0), Circle(Point2f(5,5), 0.5f0)]
        SpatialShapes(geoms, DataFrame(x_centroid=[1.0,5.0], y_centroid=[1.0,5.0], cell_id=1:2),
                      Dict{String,Any}())
    end
    function _make_tbl()
        data = [10.0 0.0; 0.0 8.0; 5.0 5.0]
        obs  = DataFrame(cell_id = 1:3)
        var  = DataFrame(gene = ["GeneA", "GeneB"])
        SpatialTable(data, obs, var, Dict{String,Any}())
    end

    # -------------------------------------------------------------------------
    @testset "crop — rectangular extent" begin
        pts = _make_pts()
        shp = _make_shp()
        ext = SpatialExtent(0.0, 3.0, 0.0, 3.0)

        s_pts = SpatialAnalysis.crop(pts, ext)
        @test size(s_pts.coordinates, 1) == 2
        @test s_pts.features.gene == ["A", "B"]

        s_shp = SpatialAnalysis.crop(shp, ext)
        @test length(s_shp.shapes) == 1
        @test Tables.getcolumn(s_shp, :cell_id) == [1]
    end

    @testset "crop — polygon ROI" begin
        pts = _make_pts()
        roi = SpatialShapes(
            [Polygon([Point2f(0,0), Point2f(3,0), Point2f(3,3), Point2f(0,3), Point2f(0,0)])],
            DataFrame(x_centroid=[1.5], y_centroid=[1.5]),
            Dict{String,Any}())
        s = SpatialAnalysis.crop(pts, roi)
        @test size(s.coordinates, 1) == 2
    end

    @testset "distances — to SpatialExtent" begin
        pts = _make_pts()
        ext = SpatialExtent(0.0, 10.0, 0.0, 10.0)
        d   = SpatialAnalysis.distances(pts, ext)
        @test d isa Vector{Float64}
        @test length(d) == 4
        @test d[1] ≈ 1.0   # (1,1) → nearest edge of [0,10]^2 is 1 away
        @test d[4] ≈ 4.0   # (6,6) → nearest edge is 4 away (10-6)
        # point outside extent
        pts2 = SpatialPoints(Float32[12 5; 5 -2], DataFrame(id=1:2), Dict{String,Any}())
        d2 = SpatialAnalysis.distances(pts2, ext)
        @test d2[1] ≈ 2.0  # (12,5) → 2 units right of xmax=10
        @test d2[2] ≈ 2.0  # (5,-2) → 2 units below ymin=0
    end

    @testset "distances — to SpatialShapes (unsigned)" begin
        shp = _make_shp()
        pts = _make_pts()
        d   = SpatialAnalysis.distances(pts, shp)
        @test d isa Vector{Float64}
        @test length(d) == 4
        @test all(d .>= 0.0)
        # point at (1,1) is inside circle at (1,1) r=0.5 → distance = 0
        @test d[1] ≈ 0.0
    end

    @testset "distances — to SpatialShapes (signed)" begin
        shp = _make_shp()
        pts = _make_pts()
        d   = SpatialAnalysis.distances(pts, shp; signed=true)
        @test d[1] < 0.0   # (1,1) is the centre of circle r=0.5 → inside → negative
        @test d[2] > 0.0   # (2,2) is outside both circles
    end

    @testset "distances — SpatialShapes to SpatialExtent" begin
        shp = _make_shp()
        ext = SpatialExtent(0.0, 10.0, 0.0, 10.0)
        d   = SpatialAnalysis.distances(shp, ext)
        @test d isa Vector{Float64}
        @test length(d) == 2
        @test all(d .>= 0.0)
    end

    @testset "classify_cells — marker-based" begin
        tbl     = _make_tbl()
        markers = Dict("TypeA" => ["GeneA"], "TypeB" => ["GeneB"])
        labels  = SpatialAnalysis.classify_cells(tbl, markers)
        @test labels isa Vector{String}
        @test length(labels) == 3
        @test labels[1] == "TypeA"
        @test labels[2] == "TypeB"
        # gene not in table → unknown
        labels2 = SpatialAnalysis.classify_cells(tbl, Dict("Z" => ["ZZZ"]); min_markers=1)
        @test all(==("unknown"), labels2)
    end

    @testset "classify_cells — signature-based dot product" begin
        tbl    = _make_tbl()
        sigs   = [1.0 0.0; 0.0 1.0]
        labels = SpatialAnalysis.classify_cells(tbl, sigs;
                     var_names=["GeneA","GeneB"], type_names=["TypeA","TypeB"])
        @test labels isa Vector{String}
        @test labels[1] == "TypeA"
        @test labels[2] == "TypeB"
    end

    @testset "classify_cells — signature-based correlation" begin
        tbl    = _make_tbl()
        sigs   = [1.0 0.0; 0.0 1.0]
        labels = SpatialAnalysis.classify_cells(tbl, sigs;
                     var_names=["GeneA","GeneB"], type_names=["TypeA","TypeB"],
                     method=:correlation)
        @test labels isa Vector{String}
        @test length(labels) == 3
    end

end
