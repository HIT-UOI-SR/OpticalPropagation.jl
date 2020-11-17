using OpticalPropagation
using Documenter

makedocs(;
    modules=[OpticalPropagation],
    authors="Yong-an Lu <miroox@outlook.com>",
    repo="https://github.com/HIT-UOI-SR/OpticalPropagation.jl/blob/{commit}{path}#L{line}",
    sitename="OpticalPropagation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HIT-UOI-SR.github.io/OpticalPropagation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "2D Light Field" => "lightfield2d.md",
        "Optical Propagation Methods" => "methods.md",
    ],
)

deploydocs(;
    repo="github.com/HIT-UOI-SR/OpticalPropagation.jl",
)
