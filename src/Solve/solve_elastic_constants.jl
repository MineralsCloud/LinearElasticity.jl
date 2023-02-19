using LinearElasticityBase: ElasticConstants
using LinearSolve: LinearProblem, solve

export LinearSystemMaker, solve_elastic_constants

struct LinearSystemMaker{X,Y,C<:SymmetryConstraint}
    x::Vector{X}
    y::Vector{Y}
    cons::C
    function LinearSystemMaker{X,Y,C}(𝐱, 𝐲, cons) where {X,Y,C}
        if length(𝐱) != length(𝐲)
            throw(DimensionMismatch("the lengths of strains and stresses must match!"))
        end
        N = minimal_npairs(cons)
        if length(𝐱) < N
            throw(ArgumentError("the number of strains/stresses must be at least $N."))
        end
        return new(𝐱, 𝐲, cons)
    end
end
LinearSystemMaker(
    𝐱::AbstractVector{X}, 𝐲::AbstractVector{Y}, cons::C=TriclinicConstraint()
) where {X,Y,C} = LinearSystemMaker{X,Y,C}(𝐱, 𝐲, cons)

function make(maker::LinearSystemMaker{<:EngineeringStrain,<:EngineeringStress})
    x, y, constraint = maker.x, maker.y, maker.cons
    𝐛 = vcat(y...)  # Length 6n vector, n = length(strains) = length(stresses)
    A = make_linear_operator(x, constraint)  # Size 6n×N matrix, N = # independent coefficients
    return LinearProblem(A, 𝐛)
end
make(maker::LinearSystemMaker{<:TensorStrain,<:TensorStress}) =
    make(LinearSystemMaker(to_voigt.(maker.x), to_voigt.(maker.y), maker.cons))
make(maker::LinearSystemMaker) = make(LinearSystemMaker(maker.y, maker.x, maker.cons))

construct(::Type{<:EngineeringStrain}, ::Type{<:EngineeringStress}) = construct_cᵢⱼ
construct(::Type{<:EngineeringStress}, ::Type{<:EngineeringStrain}) = construct_sᵢⱼ
construct(::Type{<:TensorStrain}, ::Type{<:TensorStress}) = StiffnessTensor
construct(::Type{<:TensorStress}, ::Type{<:TensorStrain}) = ComplianceTensor

function solve_elastic_constants(𝐱, 𝐲, cons=TriclinicConstraint())
    maker = LinearSystemMaker(𝐱, 𝐲, cons)
    problem = make(maker)
    solution = solve(problem)
    return construct(eltype(maker.x), eltype(maker.y))(solution)
end

minimal_npairs(::CubicConstraint) = 1
minimal_npairs(::HexagonalConstraint) = 2
minimal_npairs(::TrigonalConstraint) = 2
minimal_npairs(::TetragonalConstraint) = 2
minimal_npairs(::OrthorhombicConstraint) = 3
minimal_npairs(::MonoclinicConstraint) = 5
minimal_npairs(::TriclinicConstraint) = 6
