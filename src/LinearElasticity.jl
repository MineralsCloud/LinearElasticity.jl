module LinearElasticity

using Reexport: @reexport

@reexport using LinearElasticityBase

export CubicConstraint,
    HexagonalConstraint,
    TetragonalConstraint,
    TrigonalConstraint,
    OrthorhombicConstraint,
    MonoclinicConstraint,
    TriclinicConstraint

abstract type SymmetryConstraint end
struct TriclinicConstraint <: SymmetryConstraint end
struct MonoclinicConstraint <: SymmetryConstraint end
struct OrthorhombicConstraint <: SymmetryConstraint end
struct TetragonalConstraint <: SymmetryConstraint end
struct CubicConstraint <: SymmetryConstraint end
struct TrigonalConstraint <: SymmetryConstraint end
struct HexagonalConstraint <: SymmetryConstraint end

include("SymmetryCriteria.jl")
# include("StabilityCriteria.jl")
include("Isotropic.jl")
include("misc.jl")
include("Solve/Solve.jl")
include("ULICS.jl")

end
