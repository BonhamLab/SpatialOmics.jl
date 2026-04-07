using Documenter
using SpatialOmics
import SpatialOmicsBase
import SpatialIO
import SpatialViz

makedocs(
    sitename = "SpatialOmics.jl",
    authors  = "stx_dev contributors",
    modules  = [SpatialOmics, SpatialOmicsBase, SpatialIO, SpatialViz],
    remotes  = nothing,
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://Bonhamlab.github.io/SpatialOmics.jl",
    ),
    pages = [
        "Home"      => "index.md",
        "Guides"    => [
            "guides/quickstart.md",
            "guides/visualization.md",
            "guides/cosmx.md",
        ],
        "Reference" => [
            "reference/data_structures.md",
            "reference/views.md",
            "reference/platform_readers.md",
            "reference/spatialdata_io.md",
            "reference/format_utils.md",
        ],
    ],
    checkdocs = :exports,
    warnonly  = [:missing_docs, :docs_block, :cross_references],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo   = "github.com/Bonhamlab/SpatialOmics.jl.git",
        target = "build",
        push_preview = true,
    )
end
