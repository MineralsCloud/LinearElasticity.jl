using LinearElasticity
using Documenter

makedocs(;
    modules = [LinearElasticity],
    authors = "Qi Zhang <singularitti@outlook.com>",
    repo = "https://github.com/singularitti/LinearElasticity.jl/blob/{commit}{path}#L{line}",
    sitename = "LinearElasticity.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://singularitti.github.io/LinearElasticity.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/singularitti/LinearElasticity.jl")
