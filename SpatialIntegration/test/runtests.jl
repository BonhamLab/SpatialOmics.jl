using Test
using SpatialIntegration
using SpatialOmics

@testset "SpatialIntegration.jl" begin

    @testset "Exported symbols present" begin
        for sym in [
            :harmonize_datasets, :align_coordinate_systems, :standardize_features,
            :integrate_reference, :transfer_labels, :deconvolute_spots,
            :integrate_modalities, :joint_embedding,
        ]
            @test isdefined(SpatialIntegration, sym)
        end
    end

    @testset "Stub functions raise NotImplemented" begin
        ds1 = spatial_dataset()
        ds2 = spatial_dataset()

        @test_throws ErrorException harmonize_datasets([ds1, ds2])
        @test_throws ErrorException align_coordinate_systems([ds1, ds2])
        @test_throws ErrorException standardize_features([ds1, ds2])
        @test_throws ErrorException integrate_reference(ds1, ds2)
        @test_throws ErrorException transfer_labels(ds1, ds2)
        @test_throws ErrorException deconvolute_spots(ds1, ds2)
        @test_throws ErrorException integrate_modalities([ds1, ds2])
        @test_throws ErrorException joint_embedding([ds1, ds2])
    end

end
