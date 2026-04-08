# test/spatialdata_integration.jl
#
# Integration tests for from_spatialdata / to_spatialdata against the real
# test zarr stores in test_data/experiments/.
#
# Xenium store exercises all five element types (images, labels, points,
# shapes, tables).  The round-trip test uses a tiny synthetic dataset so it
# completes quickly regardless of host memory.

const XENIUM_ZARR = abspath(@__DIR__, "..", "..", "test_data", "experiments", "xenium_ex.zarr")
const VISIUM_ZARR = abspath(@__DIR__, "..", "..", "test_data", "experiments", "visium_ex.zarr")

# ─── helpers ─────────────────────────────────────────────────────────────────

# Check whether a value is a lazy disk-backed array (does not materialise data)
_is_lazy(x) = x isa DiskArrays.AbstractDiskArray

# ─── from_spatialdata — data-dependent tests (skipped if zarr stores absent) ──

if !isdir(XENIUM_ZARR)
    @warn "Skipping Xenium integration tests — store not found at $XENIUM_ZARR"
else

@testset "from_spatialdata — Xenium store structure" begin
    ds = from_spatialdata(XENIUM_ZARR)
    @test ds isa SpatialDataset

    # Images
    @test haskey(ds.images, "he_image")
    @test haskey(ds.images, "morphology_focus")
    @test all(img -> img isa SpatialImage, values(ds.images))

    # Labels
    @test haskey(ds.labels, "cell_labels")
    @test haskey(ds.labels, "nucleus_labels")
    @test all(lbl -> lbl isa SpatialLabels, values(ds.labels))

    # Points
    @test haskey(ds.points, "transcripts")
    @test all(pts -> pts isa SpatialPoints, values(ds.points))

    # Shapes
    @test haskey(ds.shapes, "cell_boundaries")
    @test haskey(ds.shapes, "cell_circles")
    @test haskey(ds.shapes, "nucleus_boundaries")
    @test all(shp -> shp isa SpatialShapes, values(ds.shapes))

    # Tables
    @test haskey(ds.tables, "table")
    @test all(tbl -> tbl isa SpatialTable, values(ds.tables))
end

@testset "from_spatialdata — Xenium images are lazy" begin
    ds = from_spatialdata(XENIUM_ZARR)

    img = ds.images["he_image"]
    @test _is_lazy(img.data)
    @test ndims(img.data) == 3
    @test size(img.data, 1) == 3   # c, y, x layout
    @test :x ∈ keys(img.axes)
    @test :y ∈ keys(img.axes)

    lbl = ds.labels["cell_labels"]
    @test _is_lazy(lbl.data)
    @test ndims(lbl.data) == 2    # y × x label raster
end

@testset "from_spatialdata — Xenium points dimensions" begin
    ds = from_spatialdata(XENIUM_ZARR)
    pts = ds.points["transcripts"]

    @test size(pts.coordinates, 2) == 3
    @test size(pts.coordinates, 1) == 12_165_021
    @test "feature_name" ∈ names(pts.features)
    @test "cell_id"      ∈ names(pts.features)
    @test haskey(pts.metadata, "coord_cols")
    @test pts.metadata["coord_cols"] == ["x", "y", "z"]
end

@testset "from_spatialdata — Xenium shapes" begin
    ds = from_spatialdata(XENIUM_ZARR)
    cb = ds.shapes["cell_boundaries"]
    @test length(cb.geometries) == 162_254
end

@testset "from_spatialdata — Xenium table dimensions" begin
    ds = from_spatialdata(XENIUM_ZARR)
    tbl = ds.tables["table"]

    @test size(tbl.data) == (162_254, 377)
    @test nrow(tbl.obs) == 162_254
    @test "cell_id" ∈ names(tbl.obs)
    @test nrow(tbl.var) == 377
end

end # if isdir(XENIUM_ZARR)

if !isdir(VISIUM_ZARR)
    @warn "Skipping Visium integration tests — store not found at $VISIUM_ZARR"
else

@testset "from_spatialdata — Visium store zarr layout" begin
    # The three Visium tables each hold a 5.5M × 19K sparse matrix (~1 GB each);
    # loading them in CI is impractical.  We verify the on-disk layout here and
    # leave the full from_spatialdata call to a dedicated benchmark.
    @test isfile(joinpath(VISIUM_ZARR, "zarr.json"))
    @test isdir(joinpath(VISIUM_ZARR, "images", "Visium_HD_Mouse_Small_Intestine_hires_image"))
    @test isdir(joinpath(VISIUM_ZARR, "shapes", "Visium_HD_Mouse_Small_Intestine_square_002um"))
    @test isdir(joinpath(VISIUM_ZARR, "tables", "square_002um"))
    @test isdir(joinpath(VISIUM_ZARR, "tables", "square_008um"))
    @test isdir(joinpath(VISIUM_ZARR, "tables", "square_016um"))
    # Confirm all four image groups are present (groups have zarr.json; zarr.json itself is not a group)
    img_groups = filter(e -> isdir(joinpath(VISIUM_ZARR, "images", e)), readdir(joinpath(VISIUM_ZARR, "images")))
    @test length(img_groups) == 4
end

end # if isdir(VISIUM_ZARR)

# ─── to_spatialdata / round-trip ──────────────────────────────────────────────

@testset "to_spatialdata — writes root zarr.json" begin
    mktempdir() do dir
        ds = spatial_dataset()
        to_spatialdata(ds, joinpath(dir, "empty.zarr"))
        @test isfile(joinpath(dir, "empty.zarr", "zarr.json"))
    end
end

@testset "to_spatialdata / from_spatialdata — points round-trip" begin
    mktempdir() do dir
        out = joinpath(dir, "rt.zarr")

        # Build a small dataset with 100 points in 2-D
        coords = rand(Float32, 100, 2)
        feats  = DataFrame(gene = ["g$i" for i in 1:100], cell = 1:100)
        pts    = SpatialPoints(coords, feats, Dict{String,Any}("coord_cols" => ["x", "y"]))
        ds_out = spatial_dataset()
        ds_out["spots"] = pts

        to_spatialdata(ds_out, out)
        @test isfile(joinpath(out, "points", "spots", "zarr.json"))
        @test isfile(joinpath(out, "points", "spots", "points.parquet"))

        ds_in = from_spatialdata(out)
        @test haskey(ds_in.points, "spots")

        pts_in = ds_in.points["spots"]
        @test size(pts_in.coordinates) == (100, 2)
        @test pts_in.coordinates ≈ coords
        @test "gene" ∈ names(pts_in.features)
        @test "cell" ∈ names(pts_in.features)
    end
end

@testset "to_spatialdata / from_spatialdata — image round-trip" begin
    mktempdir() do dir
        out = joinpath(dir, "rt_img.zarr")

        # 2-channel 8×8 uint8 image
        data = rand(UInt8, 2, 8, 8)
        img  = SpatialImage(data, (c = "", y = "", x = ""), Dict{String,Any}())
        ds_out = spatial_dataset()
        ds_out["dapi"] = img

        to_spatialdata(ds_out, out)
        @test isdir(joinpath(out, "images", "dapi"))

        ds_in = from_spatialdata(out)
        @test haskey(ds_in.images, "dapi")

        img_in = ds_in.images["dapi"]
        @test _is_lazy(img_in.data)
        @test size(img_in.data) == (2, 8, 8)
        @test collect(img_in.data) == data
    end
end

@testset "to_spatialdata / from_spatialdata — table round-trip" begin
    mktempdir() do dir
        out = joinpath(dir, "rt_tbl.zarr")

        X   = sparse([1, 2, 3], [1, 2, 3], Float32[1.0, 2.0, 3.0], 4, 5)
        obs = DataFrame(_index = ["c1","c2","c3","c4"], batch = ["A","A","B","B"])
        var = DataFrame(_index = ["g1","g2","g3","g4","g5"])
        tbl = SpatialTable(X, obs, var, Dict{String,Any}())
        ds_out = spatial_dataset()
        ds_out["counts"] = tbl

        to_spatialdata(ds_out, out)
        ds_in = from_spatialdata(out)
        @test haskey(ds_in.tables, "counts")

        tbl_in = ds_in.tables["counts"]
        @test size(tbl_in.data) == (4, 5)
        @test Matrix(tbl_in.data) ≈ Matrix(X)
        @test nrow(tbl_in.obs) == 4
        @test nrow(tbl_in.var) == 5
    end
end
