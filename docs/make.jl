using LinearElasticity
using Documenter

DocMeta.setdocmeta!(LinearElasticity, :DocTestSetup, :(using LinearElasticity); recursive=true)

makedocs(;
    modules=[LinearElasticity],
    authors="singularitti <singularitti@outlook.com> and contributors",
    sitename="LinearElasticity.jl",
    format=Documenter.HTML(;
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
