using Documenter, LightPropagation

makedocs(;
    modules=[LightPropagation],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/miRoox/LightPropagation.jl/blob/{commit}{path}#L{line}",
    sitename="LightPropagation.jl",
    authors="miRoox",
    assets=String[],
)

deploydocs(;
    repo="github.com/miRoox/LightPropagation.jl",
)
