using Documenter
using SpatialOmics

makedocs(
    sitename = "SpatialOmics.jl",
    authors  = "Kevin Bonham",
    modules  = [SpatialOmics],
    remotes  = nothing,
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://BonhamLab.github.io/SpatialOmics.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Introduction" => [
            "explanation/data_model.md",
            "explanation/coordinate_systems.md",
            "explanation/lazy_views.md",
        ],
        "Reference" => [
            "reference/dataset.md",
            "reference/elements.md",
            "reference/coordinate_systems.md",
            "reference/views.md",
            "reference/relations_analysis.md",
            "reference/io.md",
            "reference/visualization.md",
        ],
        "Guides" => [
            "guides/quickstart.md",
            "guides/visualization.md",
            "guides/cosmx.md",
        ],
    ],
    checkdocs = :exports,
    warnonly  = [:missing_docs, :docs_block, :cross_references],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo         = "github.com/BonhamLab/SpatialOmics.jl.git",
        target       = "build",
        push_preview = true,
    )
end
