module LinearElasticity

using Reexport: @reexport

@reexport using LinearElasticityBase

export Cubic, Hexagonal, Tetragonal, Trigonal, Orthorhombic, Monoclinic, Triclinic

abstract type SymmetryConstraint end
struct Triclinic <: SymmetryConstraint end
struct Monoclinic <: SymmetryConstraint end
struct Orthorhombic <: SymmetryConstraint end
struct Tetragonal <: SymmetryConstraint end
struct Cubic <: SymmetryConstraint end
struct Trigonal <: SymmetryConstraint end
struct Hexagonal <: SymmetryConstraint end

include("Symmetry.jl")
# include("StabilityCriteria.jl")
include("Isotropic.jl")
include("misc.jl")
include("Solve/Solve.jl")
include("ULICS.jl")

end
