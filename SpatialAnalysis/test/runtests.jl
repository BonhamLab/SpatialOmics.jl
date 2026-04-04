using Test
using SpatialAnalysis
using SpatialOmics

@testset "SpatialAnalysis.jl" begin

    @testset "Sub-module availability" begin
        # Verify sub-modules are accessible
        @test isdefined(SpatialAnalysis, :QualityControl)
        @test isdefined(SpatialAnalysis, :Preprocessing)
        @test isdefined(SpatialAnalysis, :SpatialStatistics)
        @test isdefined(SpatialAnalysis, :CellInteraction)
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
        # Statistics
        @test isdefined(SpatialAnalysis, :spatial_autocorrelation)
        @test isdefined(SpatialAnalysis, :hotspot_detection)
        @test isdefined(SpatialAnalysis, :spatial_clustering)
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

end
