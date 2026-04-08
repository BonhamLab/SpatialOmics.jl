using Test
using SpatialOmicsBase
using DiskArrays
using DataFrames
using SparseArrays
import SpatialIO
import SpatialIO:
    PlatformReader,
    XeniumReader, VisiumReader, CosMxReader, MerfishReader,
    validate_path, read_data,
    from_spatialdata, to_spatialdata,
    write_hdf5, read_hdf5, write_zarr, open_zarr

const XENIUM_MOCK = joinpath(@__DIR__, "fixtures", "xenium_mock")
const VISIUM_MOCK = joinpath(@__DIR__, "fixtures", "visium_mock")
const COSMX_MOCK  = joinpath(@__DIR__, "fixtures", "cosmx_mock")

@testset "SpatialIO.jl" begin

    @testset "PlatformReader abstract type" begin
        @test XeniumReader <: PlatformReader
        @test VisiumReader <: PlatformReader
        @test CosMxReader  <: PlatformReader
        @test MerfishReader <: PlatformReader
    end

    @testset "Reader construction (kwarg defaults)" begin
        xr = XeniumReader()
        @test xr.lazy == true
        @test xr.chunk_size == 100_000

        vr = VisiumReader()
        @test vr.lazy == true
        @test vr.image_key == "hires"

        cr = CosMxReader()
        @test cr.lazy == true
        @test isnothing(cr.fov_subset)

        mr = MerfishReader()
        @test mr.lazy == true
        @test mr.z_slice == 2
    end

    @testset "XeniumReader.validate_path with mock data" begin
        xr = XeniumReader()
        # The mock directory has the required sentinel files
        @test validate_path(xr, XENIUM_MOCK)
        # A non-existent path should fail validation
        @test !validate_path(xr, "/nonexistent/path")
    end

    @testset "VisiumReader.validate_path with mock data" begin
        vr = VisiumReader()
        @test validate_path(vr, VISIUM_MOCK)
        @test !validate_path(vr, "/nonexistent/path")
    end

    @testset "read_data raises NotImplemented (stub guard)" begin
        xr = XeniumReader()
        @test_throws ErrorException read_data(xr, XENIUM_MOCK)

        vr = VisiumReader()
        @test_throws ErrorException read_data(vr, VISIUM_MOCK)
    end

    @testset "CosMxReader.validate_path with mock data" begin
        cr = CosMxReader()
        # Accepts path containing DecodedFiles/ as direct child
        @test validate_path(cr, COSMX_MOCK)
        # Also accepts path pointing directly at DecodedFiles/
        @test validate_path(cr, joinpath(COSMX_MOCK, "DecodedFiles"))
        @test !validate_path(cr, "/nonexistent/path")
    end

    @testset "from_spatialdata — missing path returns error" begin
        @test_throws ErrorException from_spatialdata("/nonexistent/path.zarr")
    end

    @testset "to_spatialdata — writes to disk without error" begin
        mktempdir() do dir
            ds = spatial_dataset()
            @test (to_spatialdata(ds, joinpath(dir, "out.zarr")); true)
        end
    end

end

include("spatialdata_integration.jl")
include("io_backends.jl")
