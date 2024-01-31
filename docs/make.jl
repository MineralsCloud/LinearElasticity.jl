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
        "Manual" => [
            "Installation Guide" => "man/installation.md",
            "Troubleshooting" => "man/troubleshooting.md",
        ],
        "Reference" => Any[
            "Public API" => "lib/public.md",
            # "Internals" => map(
            #     s -> "lib/internals/$(s)",
            #     sort(readdir(joinpath(@__DIR__, "src/lib/internals")))
            # ),
        ],
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style-guide.md",
            "Design Principles" => "developers/design-principles.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/LinearElasticity.jl",
    devbranch="main",
)
