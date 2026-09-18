using GeometryBasics: Point2f, Polygon

@testset "Acquisition provenance" begin
    square(xmin, xmax, ymin, ymax) = Polygon([
        Point2f(xmin, ymin), Point2f(xmax, ymin), Point2f(xmax, ymax),
        Point2f(xmin, ymax), Point2f(xmin, ymin),
    ])

    function overlapping_dataset()
        ds = SpatialDataset()
        push!(ds, CoordinateSystem("global"))
        ds["fovs"] = SpatialShapes(
            [square(0, 10, 0, 10), square(5, 15, 0, 10)];
            instance_id=Int32[1, 2],
            origins=["fov_a", "fov_b"],
            coord_system="global",
        )
        push!(ds, AcquisitionSource("fov_a"; region="fovs", instance_id=1))
        push!(ds, AcquisitionSource("fov_b"; region="fovs", instance_id=2))
        ds
    end

    @testset "source membership is distinct from geometric membership" begin
        ds = overlapping_dataset()
        try
            ds["transcripts"] = SpatialPoints(
                [Point2f(7, 5), Point2f(7, 5), Point2f(12, 5)];
                origins=["fov_a", "fov_b", "fov_b"],
                coord_system="global",
            )
            ds["cells"] = SpatialShapes(
                [square(6, 8, 4, 6), square(6, 8, 4, 6)];
                origins=["fov_a", "fov_b"],
                coord_system="global",
            )

            @test sources(ds) == ["fov_a", "fov_b"]
            @test source(ds, SubString("xfov_a", 2)) == source(ds, "fov_a")
            fov_a_points = points(view(ds, "fov_a"), "transcripts")
            @test length(fov_a_points) == 1
            @test source(fov_a_points, 1) == "fov_a"
            @test origin_ids(fov_a_points) == Int32[1]
            @test length(points(view(ds, "fov_b"), "transcripts")) == 2
            @test length(shapes(view(ds, "fov_a"), "cells")) == 1
            @test length(shapes(view(ds, "fov_b"), "cells")) == 1

            overlap = SpatialExtent(6, 8, 4, 6; coord_system="global")
            @test length(points(view(ds, overlap), "transcripts")) == 2
            @test length(shapes(view(ds, overlap), "cells")) == 2
        finally
            close(ds; discard=true)
        end
    end

    @testset "missing provenance has an explicit dataset fallback" begin
        ds = overlapping_dataset()
        try
            untracked = SpatialPoints(
                [Point2f(2, 5), Point2f(12, 5)]; coord_system="global",
            )
            ds["untracked"] = untracked
            @test_logs (:warn, r"no acquisition provenance") begin
                @test length(points(view(ds, "fov_a"), "untracked")) == 1
            end
            @test_throws ArgumentError view(untracked, source(ds, "fov_a"))
        finally
            close(ds; discard=true)
        end
    end

    @testset "origin metadata survives materialisation and storage" begin
        mktempdir() do path
            ds = overlapping_dataset()
            ds["transcripts"] = SpatialPoints(
                [Point2f(2, 5), Point2f(12, 5)];
                features=(quality=Float32[0.8, 0.9],),
                origins=["fov_a", "fov_b"],
                coord_system="global",
            )
            selected = collect(points(view(ds, "fov_b"), "transcripts"))
            @test origins(selected) == ["fov_a", "fov_b"]
            @test origin_ids(selected) == Int32[2]
            @test source(selected, 1) == "fov_b"
            @test features(selected, :quality) == Float32[0.9]

            save!(ds; path)
            close(ds)
            stored = read(SpatialDataZarr(), path)
            try
                @test sources(stored) == ["fov_a", "fov_b"]
                @test source(points(stored, "transcripts"), 2) == "fov_b"
                @test length(points(view(stored, "fov_b"), "transcripts")) == 1
            finally
                close(stored)
            end
        end
    end

    @testset "source footprints crop raster elements" begin
        ds = overlapping_dataset()
        try
            ds["image"] = SpatialImage(
                reshape(Float32.(1:150), 15, 10);
                axes=(:x, :y), coord_system="global",
            )
            ds["labels"] = SpatialLabels(
                reshape(Int32.(1:150), 15, 10);
                axes=(:x, :y), coord_system="global",
            )
            @test size(images(view(ds, "fov_a"), "image")) == (10, 10)
            @test size(labels(view(ds, "fov_b"), "labels")) == (10, 10)
        finally
            close(ds; discard=true)
        end
    end

    @testset "constructor validation" begin
        @test_throws ArgumentError AcquisitionSource("bad"; region="fovs")
        @test_throws DimensionMismatch SpatialPoints(
            [Point2f(1, 1), Point2f(2, 2)]; origins=["only_one"],
        )
        @test_throws ArgumentError SpatialPoints(
            [Point2f(1, 1)]; origin_id=Int32[2], origin_codebook=["one"],
        )
    end
end
