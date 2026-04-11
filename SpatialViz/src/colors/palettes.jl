# SpatialViz/src/colors/palettes.jl
# Color palette utilities for biological spatial data visualization.
# TODO: I think all of this functionality is already defined
# in Colors.jl, ColorSchemes.jl, or Makie.jl. 
# No need to re-impliment
"""
    categorical_palette(n::Int; colorblind_safe::Bool = true) -> Vector

Return a vector of `n` distinguishable colors for categorical data.
When `colorblind_safe=true`, restricts to palettes accessible for the most
common forms of color-vision deficiency.
"""
function categorical_palette(n::Int; colorblind_safe::Bool = true)
    # TODO Phase 1: implement using ColorSchemes.jl or Makie palettes
    error("categorical_palette: not yet implemented (Phase 1 deliverable)")
end

"""
    continuous_colormap(name::Symbol = :viridis; reverse::Bool = false)

Return a continuous colormap suitable for expression or density data.
Defaults to `:viridis` which is perceptually uniform and colorblind-safe.
"""
function continuous_colormap(name::Symbol = :viridis; reverse::Bool = false)
    error("continuous_colormap: not yet implemented (Phase 1 deliverable)")
end

"""
    colorblind_safe_palette() -> Vector

Return a curated 8-color palette suitable for all common forms of
color-vision deficiency (Wong 2011 palette).
"""
function colorblind_safe_palette()
    # Wong (2011) Nature Methods colorblind-safe 8-color palette (RGB hex)
    return [
        "#000000",  # Black
        "#E69F00",  # Orange
        "#56B4E9",  # Sky blue
        "#009E73",  # Bluish green
        "#F0E442",  # Yellow
        "#0072B2",  # Blue
        "#D55E00",  # Vermillion
        "#CC79A7",  # Reddish purple
    ]
end
