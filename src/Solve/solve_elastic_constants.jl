export Problem, solve, solve_elastic_constants

struct Problem{X,Y,C<:SymmetryConstraint}
    x::Vector{X}
    y::Vector{Y}
    cons::C
    function Problem{X,Y,C}(𝐱, 𝐲, cons) where {X,Y,C}
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
Problem(
    𝐱::AbstractVector{X}, 𝐲::AbstractVector{Y}, cons::C=TriclinicConstraint()
) where {X,Y,C} = Problem{X,Y,C}(𝐱, 𝐲, cons)

function solve(problem::Problem)
    x, y, constraint = problem.x, problem.y, problem.cons
    𝐛 = vcat(y...)  # Length 6n vector, n = length(strains) = length(stresses)
    A = construct_linear(x, constraint)  # Size 6n×N matrix, N = # independent coefficients
    𝐱 = A \ 𝐛  # Length N vector
    return construct(eltype(x), eltype(y))(𝐱, constraint)
end
function solve(problem::Problem)
    constants = solve(Problem(to_voigt.(problem.x), to_voigt.(problem.y), problem.cons))
    return construct(eltype(problem.x), eltype(problem.y))(constants)
end

construct(::Type{<:EngineeringStrain}, ::Type{<:EngineeringStress}) = construct_cᵢⱼ
construct(::Type{<:EngineeringStress}, ::Type{<:EngineeringStrain}) = construct_sᵢⱼ
construct(::Type{<:TensorStrain}, ::Type{<:TensorStress}) = StiffnessTensor
construct(::Type{<:TensorStress}, ::Type{<:TensorStrain}) = ComplianceTensor

solve_elastic_constants(𝐱, 𝐲, cons=TriclinicConstraint()) = solve(Problem(𝐱, 𝐲, cons))

minimal_npairs(::CubicConstraint) = 1
minimal_npairs(::HexagonalConstraint) = 2
minimal_npairs(::TrigonalConstraint) = 2
minimal_npairs(::TetragonalConstraint) = 2
minimal_npairs(::OrthorhombicConstraint) = 3
minimal_npairs(::MonoclinicConstraint) = 5
minimal_npairs(::TriclinicConstraint) = 6
