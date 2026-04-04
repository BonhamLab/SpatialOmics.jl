using Test
using SpatialViz
using SpatialOmics

@testset "SpatialViz.jl" begin

    @testset "Backend type hierarchy" begin
        @test MakieBackend <: VisualizationBackend
        @test WGLMakieBackend <: VisualizationBackend
        @test NapariBackend <: VisualizationBackend
        # BonitoBackend is an alias for WGLMakieBackend
        @test BonitoBackend === WGLMakieBackend
    end

    @testset "Exported symbols present" begin
        for sym in [
            :spatial_plot, :plot_points, :plot_image, :plot_labels, :plot_shapes,
            :overlay_points_on_image,
            :annotation_tool, :roi_select, :lasso_select, :rectangle_select,
            :save_annotations, :load_annotations,
            :categorical_palette, :continuous_colormap, :colorblind_safe_palette,
        ]
            @test isdefined(SpatialViz, sym)
        end
    end

    @testset "colorblind_safe_palette returns 8 colors" begin
        pal = colorblind_safe_palette()
        @test length(pal) == 8
        @test all(c -> startswith(c, "#"), pal)
    end

    @testset "Stub functions raise NotImplemented" begin
        ds = spatial_dataset()
        @test_throws ErrorException spatial_plot(ds)
        @test_throws ErrorException categorical_palette(5)
        @test_throws ErrorException continuous_colormap(:viridis)
    end

end
