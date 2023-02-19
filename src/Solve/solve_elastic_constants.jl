export Problem, solve

struct Problem{X,Y,C<:SymmetryConstraint}
    x::Vector{X}
    y::Vector{Y}
    cons::C
end
Problem(𝐱, 𝐲, cons=TriclinicConstraint()) = Problem(𝐱, 𝐲, cons)

function solve(problem::Problem{<:EngineeringStress,<:EngineeringStrain})
    strains, stresses, constraint = problem.x, problem.y, problem.cons
    if length(strains) != length(stresses)
        throw(DimensionMismatch("the lengths of strains and stresses must match!"))
    end
    n = minimal_npairs(constraint)
    if length(strains) < n
        throw(ArgumentError("the number of strains/stresses must be at least $n."))
    end
    σ = vcat(stresses...)  # Length 6n vector, n = length(strains) = length(stresses)
    ε = combine_strains(strains, constraint)  # Size 6n×N matrix, N = # independent coefficients
    𝐜 = ε \ σ  # Length N vector
    return construct_cᵢⱼ(𝐜, constraint)
end
function solve(problem::Problem{<:EngineeringStrain,<:EngineeringStress})
    stresses, strains, constraint = problem.x, problem.y, problem.cons
    if length(strains) != length(stresses)
        throw(DimensionMismatch("the lengths of strains and stresses must match!"))
    end
    n = minimal_npairs(constraint)
    if length(strains) < n
        throw(ArgumentError("the number of strains/stresses must be at least $n."))
    end
    ε = vcat(strains...)
    σ = combine_stresses(stresses, constraint)
    𝐬 = σ \ ε
    return construct_sᵢⱼ(𝐬, constraint)
end
function solve(problem::Problem{<:TensorStrain,<:TensorStress})
    cᵢⱼ = solve(Problem(to_voigt.(problem.x), to_voigt.(problem.y), problem.cons))
    return StiffnessTensor(cᵢⱼ)
end
function solve(problem::Problem{<:TensorStress,<:TensorStrain})
    sᵢⱼ = solve(Problem(to_voigt.(problem.x), to_voigt.(problem.y), problem.cons))
    return ComplianceTensor(sᵢⱼ)
end

minimal_npairs(::CubicConstraint) = 1
minimal_npairs(::HexagonalConstraint) = 2
minimal_npairs(::TrigonalConstraint) = 2
minimal_npairs(::TetragonalConstraint) = 2
minimal_npairs(::OrthorhombicConstraint) = 3
minimal_npairs(::MonoclinicConstraint) = 5
minimal_npairs(::TriclinicConstraint) = 6
