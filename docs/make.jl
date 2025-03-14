using SpatialOmics
using Documenter

DocMeta.setdocmeta!(SpatialOmics, :DocTestSetup, :(using SpatialOmics); recursive=true)

makedocs(;
    modules=[SpatialOmics],
    authors="Kevin Bonham <kevin@bonham.ch> and contributors",
    sitename="SpatialOmics.jl",
    format=Documenter.HTML(;
        canonical="https://github.com/BonhamLab/SpatialOmics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kescobo/SpatialOmics.jl",
    devbranch="main",
)
