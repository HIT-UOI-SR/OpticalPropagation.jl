using Documenter, OpticalPropagation

makedocs(;
    modules=[OpticalPropagation],
    format=Documenter.HTML(assets=String[]),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/HIT-UOI-SR/OpticalPropagation.jl/blob/{commit}{path}#L{line}",
    sitename="OpticalPropagation.jl",
    authors="miRoox",
)

deploydocs(;
    repo="github.com/HIT-UOI-SR/OpticalPropagation.jl",
)
