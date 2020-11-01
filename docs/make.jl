using Documenter, OpticalPropagation

makedocs(;
    modules=[OpticalPropagation],
    format=Documenter.HTML(assets=String[], prettyurls = get(ENV, "CI", nothing) == "true"),
    pages=[
        "Home" => "index.md",
        "2D Light Field" => "lightfield2d.md",
        "Optical Propagation Methods" => "methods.md",
    ],
    repo="https://github.com/HIT-UOI-SR/OpticalPropagation.jl/blob/{commit}{path}#L{line}",
    sitename="OpticalPropagation.jl",
    authors="miRoox",
)

deploydocs(;
    repo="github.com/HIT-UOI-SR/OpticalPropagation.jl",
)
