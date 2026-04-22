using Test
using SpatialOmicsBase
using DiskArrays
using DataFrames
using SparseArrays
using GeometryBasics
import SpatialIO
import SpatialIO:
    PlatformReader,
    XeniumReader, VisiumReader, CosMxReader, MerfishReader,
    validate_path, read_data, load,
    from_spatialdata, to_spatialdata,
    write_hdf5, read_hdf5, write_zarr, open_zarr

const XENIUM_MOCK    = joinpath(@__DIR__, "fixtures", "xenium_mock")
const VISIUM_MOCK    = joinpath(@__DIR__, "fixtures", "visium_mock")
const COSMX_MOCK  = joinpath(@__DIR__, "fixtures", "cosmx_mock")

@testset "SpatialIO.jl" begin

    @testset "PlatformReader abstract type" begin
        @test XeniumReader   <: PlatformReader
        @test VisiumReader   <: PlatformReader
        @test CosMxReader <: PlatformReader
        @test MerfishReader  <: PlatformReader
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
        @test isnothing(cr.sample)
        @test isnothing(cr.morphology_dir)
        @test isnothing(cr.morphology_zarr)

        mr = MerfishReader()
        @test mr.lazy == true
        @test mr.z_slice == 2
    end

    @testset "XeniumReader.validate_path" begin
        @test  validate_path(XeniumReader(), XENIUM_MOCK)
        @test !validate_path(XeniumReader(), "/nonexistent/path")
    end

    @testset "VisiumReader.validate_path" begin
        @test  validate_path(VisiumReader(), VISIUM_MOCK)
        @test !validate_path(VisiumReader(), "/nonexistent/path")
    end

    @testset "read_data raises NotImplemented (stub guard)" begin
        @test_throws ErrorException read_data(XeniumReader(), XENIUM_MOCK)
        @test_throws ErrorException read_data(VisiumReader(), VISIUM_MOCK)
    end

    # ── CosMxReader ────────────────────────────────────────────────────────

    @testset "CosMxReader validate_path" begin
        # Root dir containing flatFiles/
        @test  validate_path(CosMxReader(), COSMX_MOCK)
        # Direct sample dir (contains *_tx_file.csv.gz)
        sample_dir = joinpath(COSMX_MOCK, "flatFiles", "mock_sample")
        @test  validate_path(CosMxReader(), sample_dir)
        @test !validate_path(CosMxReader(), "/nonexistent/path")
    end

    @testset "CosMxReader read_data" begin
        ds = read_data(CosMxReader(), COSMX_MOCK)

        # Transcripts — 3 rows in mock tx_file
        @test haskey(ds.points, "transcripts")
        pts = ds.points["transcripts"]
        @test size(pts.coordinates, 1) == 3
        @test size(pts.coordinates, 2) == 2

        # Cell boundaries — 2 cells in mock polygons file
        @test haskey(ds.shapes, "cell_boundaries")
        shp = ds.shapes["cell_boundaries"]
        @test length(shp.shapes) == 2
        # Geometries are decoded Polygons
        @test geometry(shp, 1) isa GeometryBasics.Polygon

        # Expression table — 2 cells × 3 genes
        @test haskey(ds.tables, "expression")
        tbl = ds.tables["expression"]
        @test size(tbl.data, 1) == 2
        @test size(tbl.data, 2) == 3
        @test tbl.var.gene == ["Gapdh", "Actb", "Mki67"]

        # Metadata
        @test ds.metadata["format"] == "CosMx"
        @test ds.metadata["sample"] == "mock_sample"
    end

    @testset "CosMxReader load dispatch" begin
        ds = load(CosMxReader(), COSMX_MOCK)
        @test ds isa SpatialDataset
    end

    @testset "CosMxReader explicit sample name" begin
        ds = read_data(CosMxReader(; sample="mock_sample"), COSMX_MOCK)
        @test ds isa SpatialDataset
    end

    @testset "CosMxReader bad sample name" begin
        @test_throws ErrorException read_data(CosMxReader(; sample="no_such"), COSMX_MOCK)
    end

    # ── SpatialData round-trip ────────────────────────────────────────────────

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
