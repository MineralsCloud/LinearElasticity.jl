#! format: off
using LinearElasticity
using Documenter

DocMeta.setdocmeta!(LinearElasticity, :DocTestSetup, :(using LinearElasticity); recursive=true)

makedocs(;
    modules=[LinearElasticity],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/MineralsCloud/LinearElasticity.jl/blob/{commit}{path}#{line}",
    sitename="LinearElasticity.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/LinearElasticity.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/LinearElasticity.jl",
    devbranch="main",
)
