```@meta
CurrentModule = LinearElasticity
```

# LinearElasticity

Documentation for [LinearElasticity](https://github.com/MineralsCloud/LinearElasticity.jl).

See the [Index](@ref main-index) for the complete list of documented functions
and types.

The code, which is [hosted on GitHub](https://github.com/MineralsCloud/LinearElasticity.jl), is tested
using various continuous integration services for its validity.

This repository is created and maintained by
[@singularitti](https://github.com/singularitti), and contributions are highly welcome.

## Package features

Compute and convert linear elasticity values based on the Cauchy elasticity theory.

## Installation

The package can be installed with the Julia package manager.
From [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), type `]` to enter
the [Pkg mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode) and run:

```julia-repl
pkg> add LinearElasticity
```

Or, equivalently, via [`Pkg.jl`](https://pkgdocs.julialang.org/v1/):

```@repl
import Pkg; Pkg.add("LinearElasticity")
```

## Documentation

- [**STABLE**](https://MineralsCloud.github.io/LinearElasticity.jl/stable) — **documentation of the most recently tagged version.**
- [**DEV**](https://MineralsCloud.github.io/LinearElasticity.jl/dev) — _documentation of the in-development version._

## Project status

The package is developed for and tested against Julia `v1.6` and above on Linux, macOS, and
Windows.

## Questions and contributions

You can post usage questions on
[our discussion page](https://github.com/MineralsCloud/LinearElasticity.jl/discussions).

We welcome contributions, feature requests, and suggestions. If you encounter any problems,
please open an [issue](https://github.com/MineralsCloud/LinearElasticity.jl/issues).
The [Contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

## Manual outline

```@contents
Pages = [
    "man/installation.md",
    "man/troubleshooting.md",
    "developers/contributing.md",
    "developers/style-guide.md",
    "developers/design-principles.md",
]
Depth = 3
```

## Library outline

```@contents
Pages = ["lib/public.md", "lib/internals.md"]
```

### [Index](@id main-index)

```@index
Pages = ["lib/public.md"]
```
