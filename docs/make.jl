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
        "Manual" => [
            "Installation Guide" => "installation.md",
        ],
        # "Public API" => "public.md",
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style-guide.md",
            "Design Principles" => "developers/design-principles.md",
        ],
        "Troubleshooting" => "troubleshooting.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/LinearElasticity.jl",
    devbranch="main",
)
