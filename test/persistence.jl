using GeometryBasics: Point2f

@testset "Explicit persistence" begin
    @testset "changes are staged and selectively saved" begin
        mktempdir() do path
            ds = SpatialDataset(; path)
            @test !isdirty(ds)

            push!(ds, CoordinateSystem("global"))
            ds["transcripts"] = SpatialPoints(
                [Point2f(1, 2), Point2f(3, 4)];
                features=(fov=Int32[1, 2], compartment=["nucleus", "cytoplasm"]),
                coord_system="global",
            )
            ds.metadata["sample"] = "A"

            @test isdirty(ds)
            @test Set(c.kind for c in dirty(ds)) == Set((:dataset, :element, :metadata))
            @test !isdir(joinpath(path, "points", "transcripts"))

            save!(ds, "transcripts")
            @test isdir(joinpath(path, "points", "transcripts"))
            @test isdirty(ds)
            @test all(c.name != "transcripts" for c in dirty(ds))

            save!(ds)
            @test !isdirty(ds)

            stored = read(SpatialDataZarr(), path)
            @test length(points(stored, "transcripts")) == 2
            @test features(points(stored, "transcripts"), :fov) == Int32[1, 2]
            @test features(points(stored, "transcripts"), :compartment) == ["nucleus", "cytoplasm"]
            @test "global" in coord_systems(stored)
            @test stored.metadata["sample"] == "A"
            close(stored)
            close(ds)
        end
    end

    @testset "close requires an explicit choice" begin
        ds = SpatialDataset()
        path = ds.backing.path
        ds["points"] = SpatialPoints([Point2f(1, 1)])
        @test_throws ArgumentError close(ds)
        @test isdir(path)
        close(ds; discard=true)
        @test !isdir(path)
    end

    @testset "discard restores saved state" begin
        mktempdir() do path
            ds = SpatialDataset(; path)
            ds["points"] = SpatialPoints([Point2f(1, 1), Point2f(2, 2)])
            save!(ds)

            edit!(ds, "points") do points
                points.coords[1] = Point2f(9, 9)
            end
            abandoned = points(ds, "points")
            @test isdirty(ds)
            discard!(ds, "points")
            @test !isdirty(ds)
            @test coords(points(ds, "points"))[1] == Point2f(1, 1)
            @test SpatialOmics._owning_dataset(abandoned) === nothing
            close(ds)
        end
    end

    @testset "supported mutators mark attached elements dirty" begin
        ds = SpatialDataset()
        pts = SpatialPoints([Point2f(1, 1)]; coord_system="local")
        ds["points"] = pts
        save!(ds)
        apply!(SpatialOmics.translation(1, 2, "local", "global"), pts)
        @test isdirty(ds)
        @test coord_system(pts) == "global"
        close(ds; discard=true)
    end

    @testset "sequence transformations round-trip" begin
        mktempdir() do path
            ds = SpatialDataset(; path)
            push!(ds, CoordinateSystem("a"))
            push!(ds, CoordinateSystem("b"))
            push!(ds, CoordinateSystem("c"))
            first_step = SpatialOmics.translation(1, 0, "a", "b")
            second_step = SpatialOmics.translation(0, 2, "b", "c")
            push!(ds, Sequence([first_step, second_step], "a", "c"))
            save!(ds)

            stored = read(SpatialDataZarr(), path)
            @test only(stored.transforms) isa Sequence
            @test apply(only(stored.transforms), [0.0 0.0]) ≈ [1.0 2.0]
            close(stored)
            close(ds)
        end
    end

    @testset "deletion and cross-type replacement remove stale storage" begin
        mktempdir() do path
            ds = SpatialDataset(; path)
            ds["object"] = SpatialPoints([Point2f(1, 1)])
            save!(ds)
            @test isdir(joinpath(path, "points", "object"))

            ds["object"] = SpatialShapes(
                [Polygon([Point2f(0, 0), Point2f(2, 0), Point2f(2, 2), Point2f(0, 0)])],
            )
            save!(ds)
            @test !isdir(joinpath(path, "points", "object"))
            @test isdir(joinpath(path, "shapes", "object"))

            delete!(ds, "object")
            @test any(c.state == :deleted for c in dirty(ds))
            save!(ds)
            @test !isdir(joinpath(path, "shapes", "object"))
            close(ds)
        end
    end

    @testset "write! remains a saving compatibility API" begin
        mktempdir() do parent
            ds = SpatialDataset()
            old_path = ds.backing.path
            ds["points"] = SpatialPoints([Point2f(1, 1)])
            target = joinpath(parent, "dataset.zarr")
            write!(ds, target, SpatialDataZarr())

            @test ds.backing.path == abspath(target)
            @test !ds.backing.owned
            @test !isdirty(ds)
            @test !isdir(old_path)
            @test isdir(joinpath(target, "points", "points"))
            close(ds)
        end
    end

    @testset "temporary stores reject nested save targets" begin
        ds = SpatialDataset()
        ds["points"] = SpatialPoints([Point2f(1, 1)])
        nested = joinpath(ds.backing.path, "nested.zarr")
        @test_throws ArgumentError save!(ds; path=nested)
        @test isdirty(ds)
        close(ds; discard=true)
    end

    @testset "unsupported attachments fail immediately" begin
        ds = SpatialDataset()
        @test_throws MethodError setindex!(ds, (x=1,), "unsupported")
        close(ds)
    end

    @testset "unsupported metadata remains dirty after a failed save" begin
        ds = SpatialDataset()
        ds.metadata["bad"] = (objects=Any[Ref(1)],)
        @test_throws ArgumentError save!(ds)
        @test isdirty(ds)
        close(ds; discard=true)
    end

    @testset "display reports unsaved state" begin
        ds = SpatialDataset()
        ds["points"] = SpatialPoints([Point2f(1, 1)])
        @test occursin("unsaved", sprint(show, ds))
        @test occursin("unsaved", sprint(show, MIME("text/plain"), ds))
        close(ds; discard=true)
    end

    @testset "native stores require an explicit format version" begin
        mktempdir() do path
            open(joinpath(path, "zarr.json"), "w") do io
                write(io, """{"zarr_format":3,"node_type":"group","attributes":{}}""")
            end
            open(joinpath(path, "spatialomics_meta.json"), "w") do io
                write(io, """{"coord_systems":[]}""")
            end
            @test native_store_version(path) === nothing
            error = try
                read(SpatialDataZarr(), path)
                nothing
            catch exception
                exception
            end
            @test error isa ArgumentError
            @test occursin("rebuild", sprint(showerror, error))
            @test occursin("does not upgrade stores automatically", sprint(showerror, error))
        end
    end
end
