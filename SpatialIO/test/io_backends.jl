# test/io_backends.jl
#
# Tests for read_hdf5 / write_hdf5 and open_zarr / write_zarr.

using Test
using SpatialIO
using SpatialOmicsBase
using DataFrames
using SparseArrays
using Tables

# ─── helpers ─────────────────────────────────────────────────────────────────

function _make_test_dataset()
    ds = spatial_dataset()

    # image
    img_data = rand(Float32, 3, 16, 16)
    img = SpatialImage(img_data, (c = "", y = "µm", x = "µm"), Dict{String,Any}())
    ds["morphology"] = img

    # labels
    lbl_data = Int32.(rand(1:5, 16, 16))
    ds["cell_labels"] = SpatialLabels(lbl_data, Dict{String,Any}())

    # points
    coords = rand(Float32, 20, 2)
    feats  = DataFrame(gene = fill("GAPDH", 20), cell_id = Int32.(1:20))
    pts    = SpatialPoints(coords, feats,
                           Dict{String,Any}("coord_cols" => ["x", "y"]))
    ds["transcripts"] = pts

    # shapes (no geometry needed for HDF5 round-trip; nothing geometry placeholder)
    shp_feats = DataFrame(cell_id = Int32.(1:5), area = rand(Float32, 5))
    ds["cell_boundaries"] = SpatialShapes(fill(nothing, 5), shp_feats, Dict{String,Any}())

    # table
    X   = sparse(rand(Float32, 10, 50))
    obs = DataFrame(barcode = ["BC$(i)" for i in 1:10],
                    n_genes  = Int32.(rand(100:500, 10)))
    var = DataFrame(gene = ["G$(i)" for i in 1:50])
    ds["expression"] = SpatialTable(X, obs, var, Dict{String,Any}())

    return ds
end

# ─── HDF5 backend ─────────────────────────────────────────────────────────────

@testset "write_hdf5 / read_hdf5 — round-trip" begin
    mktempdir() do dir
        path = joinpath(dir, "test.h5")
        ds   = _make_test_dataset()

        write_hdf5(ds, path)
        @test isfile(path)

        ds2 = read_hdf5(path)
        @test ds2 isa SpatialDataset

        @testset "images" begin
            @test haskey(ds2.images, "morphology")
            img_orig = ds.images["morphology"]
            img_rt   = ds2.images["morphology"]
            @test size(img_rt.data) == size(img_orig.data)
            @test collect(img_rt.data) ≈ collect(img_orig.data)
            @test keys(img_rt.axes) == keys(img_orig.axes)
        end

        @testset "labels" begin
            @test haskey(ds2.labels, "cell_labels")
            @test collect(ds2.labels["cell_labels"].data) ==
                  collect(ds.labels["cell_labels"].data)
        end

        @testset "points" begin
            @test haskey(ds2.points, "transcripts")
            pts_rt = ds2.points["transcripts"]
            @test size(pts_rt.coordinates) == size(ds.points["transcripts"].coordinates)
            @test pts_rt.coordinates ≈ ds.points["transcripts"].coordinates
            @test nrow(pts_rt.features) == 20
            @test pts_rt.features.gene == ds.points["transcripts"].features.gene
        end

        @testset "shapes" begin
            @test haskey(ds2.shapes, "cell_boundaries")
            shp_rt = ds2.shapes["cell_boundaries"]
            @test length(shp_rt.shapes) == 5
            @test Tables.getcolumn(shp_rt, :cell_id) ==
                  Tables.getcolumn(ds.shapes["cell_boundaries"], :cell_id)
        end

        @testset "tables" begin
            @test haskey(ds2.tables, "expression")
            tbl_rt = ds2.tables["expression"]
            @test size(tbl_rt.data) == (10, 50)
            @test tbl_rt.obs.barcode == ds.tables["expression"].obs.barcode
            @test tbl_rt.var.gene    == ds.tables["expression"].var.gene
        end
    end
end

@testset "read_hdf5 — missing file throws" begin
    @test_throws ErrorException read_hdf5("/nonexistent/path.h5")
end

# ─── Zarr backend (wrappers) ──────────────────────────────────────────────────

@testset "open_zarr / write_zarr — round-trip (points)" begin
    mktempdir() do dir
        path = joinpath(dir, "test.zarr")
        ds   = spatial_dataset()
        coords = rand(Float32, 15, 2)
        feats  = DataFrame(gene = fill("ACTB", 15))
        ds["tx"] = SpatialPoints(coords, feats,
                                Dict{String,Any}("coord_cols" => ["x","y"]))

        write_zarr(ds, path)
        @test isdir(path)

        ds2 = open_zarr(path)
        @test ds2 isa SpatialDataset
        @test haskey(ds2.points, "tx")
        @test size(ds2.points["tx"].coordinates) == (15, 2)
    end
end
