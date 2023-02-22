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

export LinearSystemMaker, make, solve_elastic_constants

struct LinearSystemMaker{X,Y,C<:SymmetryConstraint}
    x::Vector{X}
    y::Vector{Y}
    cstr::C
    function LinearSystemMaker{X,Y,C}(𝐱, 𝐲, cstr) where {X,Y,C}
        if length(𝐱) != length(𝐲)
            throw(DimensionMismatch("the lengths of strains and stresses must match!"))
        end
        N = minimal_npairs(cstr)
        if length(𝐱) < N
            throw(ArgumentError("the number of strains/stresses must be at least $N."))
        end
        return new(𝐱, 𝐲, cstr)
    end
end
LinearSystemMaker(
    𝐱::AbstractVector{X}, 𝐲::AbstractVector{Y}, cstr::C=Triclinic()
) where {X,Y,C} = LinearSystemMaker{X,Y,C}(𝐱, 𝐲, cstr)

function make(maker::LinearSystemMaker{<:EngineeringStrain,<:EngineeringStress})
    x, y, cstr = maker.x, maker.y, maker.cstr
    𝐛 = mapreduce(collect, vcat, y)  # Length 6n vector, n = length(strains) = length(stresses)
    A = make_linear_operator(x, cstr)  # Size 6n×N matrix, N = # independent coefficients
    return LinearProblem(A, 𝐛)
end
make(maker::LinearSystemMaker{<:TensorStrain,<:TensorStress}) =
    make(LinearSystemMaker(to_voigt.(maker.x), to_voigt.(maker.y), maker.cstr))
make(maker::LinearSystemMaker) = make(LinearSystemMaker(maker.y, maker.x, maker.cstr))

target(maker::LinearSystemMaker{<:EngineeringStrain,<:EngineeringStress}) =
    Base.Fix2(construct_cᵢⱼ, maker.cstr)
target(maker::LinearSystemMaker{<:EngineeringStress,<:EngineeringStrain}) =
    Base.Fix2(construct_sᵢⱼ, maker.cstr)
target(maker::LinearSystemMaker{<:TensorStrain,<:TensorStress}) =
    StiffnessTensor ∘ Base.Fix2(construct_cᵢⱼ, maker.cstr)
target(maker::LinearSystemMaker{<:TensorStress,<:TensorStrain}) =
    ComplianceTensor ∘ Base.Fix2(construct_cᵢⱼ, maker.cstr)

function solve_elastic_constants(𝐱, 𝐲, cons=Triclinic(), args...; kwargs...)
    maker = LinearSystemMaker(𝐱, 𝐲, cons)
    problem = make(maker)
    solution = solve(problem, args...; kwargs...)
    return target(maker)(solution)
end

minimal_npairs(::Cubic) = 1
minimal_npairs(::Hexagonal) = 2
minimal_npairs(::Trigonal) = 2
minimal_npairs(::Tetragonal) = 2
minimal_npairs(::Orthorhombic) = 3
minimal_npairs(::Monoclinic) = 5
minimal_npairs(::Triclinic) = 6
