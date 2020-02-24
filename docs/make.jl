using Elasticity
using Documenter

makedocs(;
    modules=[Elasticity],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/singularitti/Elasticity.jl/blob/{commit}{path}#L{line}",
    sitename="Elasticity.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singularitti.github.io/Elasticity.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/singularitti/Elasticity.jl",
)
