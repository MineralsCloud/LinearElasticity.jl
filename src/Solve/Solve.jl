module Solve

using LinearElasticityBase:
    ElasticConstants,
    EngineeringStrain,
    EngineeringStress,
    TensorStrain,
    TensorStress,
    StiffnessMatrix,
    ComplianceMatrix,
    StiffnessTensor,
    ComplianceTensor,
    to_voigt
using LinearSolve: LinearProblem, solve

using ..Symmetry:
    SymmetryConstraint,
    Cubic,
    Hexagonal,
    Tetragonal,
    Trigonal,
    Orthorhombic,
    Monoclinic,
    Triclinic

export ProblemMaker, make, solve_elastic_constants

struct ProblemMaker{X,Y,C<:SymmetryConstraint}
    x::Vector{X}
    y::Vector{Y}
    cstr::C
    function ProblemMaker{X,Y,C}(𝐱, 𝐲, cstr) where {X,Y,C}
        if length(𝐱) != length(𝐲)
            throw(DimensionMismatch("the lengths of strains and stresses must match!"))
        end
        n = minimal_npairs(cstr)
        if length(𝐱) < n
            throw(ArgumentError("the number of strains/stresses must be at least $n."))
        end
        return new(𝐱, 𝐲, cstr)
    end
end
ProblemMaker(
    𝐱::AbstractVector{X}, 𝐲::AbstractVector{Y}, cstr::C=Triclinic()
) where {X,Y,C} = ProblemMaker{X,Y,C}(𝐱, 𝐲, cstr)

function make(maker::ProblemMaker{<:EngineeringStrain,<:EngineeringStress})
    𝐱, 𝐲, cstr = maker.x, maker.y, maker.cstr
    𝐛 = mapreduce(collect, vcat, 𝐲)  # Length 6n vector, n = length(strains) = length(stresses)
    A = make_linear_operator(𝐱, cstr)  # Size 6n×N matrix, N = # independent coefficients
    return LinearProblem(A, 𝐛)  # Linear system A 𝐮 = 𝐛
end
make(maker::ProblemMaker{<:TensorStrain,<:TensorStress}) =
    make(ProblemMaker(to_voigt.(maker.x), to_voigt.(maker.y), maker.cstr))
make(maker::ProblemMaker) = make(ProblemMaker(maker.y, maker.x, maker.cstr))

function solve_elastic_constants(𝐱, 𝐲, cstr=Triclinic(), args...; kwargs...)
    maker = ProblemMaker(𝐱, 𝐲, cstr)
    problem = make(maker)
    solution = solve(problem, args...; kwargs...)
    return target(maker)(solution)
end

target(maker::ProblemMaker{<:EngineeringStrain,<:EngineeringStress}) =
    Base.Fix2(construct_cᵢⱼ, maker.cstr)
target(maker::ProblemMaker{<:EngineeringStress,<:EngineeringStrain}) =
    Base.Fix2(construct_sᵢⱼ, maker.cstr)
target(maker::ProblemMaker{<:TensorStrain,<:TensorStress}) =
    StiffnessTensor ∘ Base.Fix2(construct_cᵢⱼ, maker.cstr)
target(maker::ProblemMaker{<:TensorStress,<:TensorStrain}) =
    ComplianceTensor ∘ Base.Fix2(construct_sᵢⱼ, maker.cstr)

minimal_npairs(::Cubic) = 1
minimal_npairs(::Hexagonal) = 2
minimal_npairs(::Trigonal) = 2
minimal_npairs(::Tetragonal) = 2
minimal_npairs(::Orthorhombic) = 3
minimal_npairs(::Monoclinic) = 5
minimal_npairs(::Triclinic) = 6

include("utils.jl")

end
