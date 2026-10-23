using OpticalPropagation
using Documenter

makedocs(;
    modules=[OpticalPropagation],
    authors="Yong-an Lu <miroox@outlook.com>",
    repo=Documenter.Remotes.GitHub("HIT-UOI-SR", "OpticalPropagation.jl"),
    sitename="OpticalPropagation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        canonical="https://HIT-UOI-SR.github.io/OpticalPropagation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "2D Light Field" => "lightfield2d.md",
        "Optical Propagation Methods" => "methods.md",
    ],
)
