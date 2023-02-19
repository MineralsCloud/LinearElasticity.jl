export solve_elastic_constants

function solve_elastic_constants(
    system::CrystalSystem,
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    if length(strains) != length(stresses)
        throw(DimensionMismatch("the lengths of strains and stresses must match!"))
    end
    n = minimal_npairs(system)
    if length(strains) < n
        throw(ArgumentError("the number of strains/stresses must be at least $n."))
    end
    σ = vcat(stresses...)  # Length 6n vector, n = length(strains) = length(stresses)
    ε = combine_strains(system, strains)  # Size 6n×N matrix, N = # independent coefficients
    𝐜 = ε \ σ  # Length N vector
    return construct_cᵢⱼ(system, 𝐜)
end
function solve_elastic_constants(
    system::CrystalSystem,
    stresses::AbstractVector{<:EngineeringStress},
    strains::AbstractVector{<:EngineeringStrain},
)
    if length(strains) != length(stresses)
        throw(DimensionMismatch("the lengths of strains and stresses must match!"))
    end
    n = minimal_npairs(system)
    if length(strains) < n
        throw(ArgumentError("the number of strains/stresses must be at least $n."))
    end
    ε = vcat(strains...)
    σ = combine_stresses(system, stresses)
    𝐬 = σ \ ε
    return construct_sᵢⱼ(system, 𝐬)
end
function solve_elastic_constants(
    system::CrystalSystem,
    strains::AbstractVector{<:TensorStrain},
    stresses::AbstractVector{<:TensorStress},
)
    cᵢⱼ = solve_elastic_constants(
        system, EngineeringStrain.(strains), EngineeringStress.(stresses)
    )
    return StiffnessTensor(cᵢⱼ)
end
function solve_elastic_constants(
    system::CrystalSystem,
    stresses::AbstractVector{<:TensorStress},
    strains::AbstractVector{<:TensorStrain},
)
    sᵢⱼ = solve_elastic_constants(
        system, EngineeringStrain.(strains), EngineeringStress.(stresses)
    )
    return ComplianceTensor(sᵢⱼ)
end
solve_elastic_constants(strains_or_stresses, stresses_or_strains) =
    solve_elastic_constants(TriclinicConstraint(), strains_or_stresses, stresses_or_strains)

minimal_npairs(::CubicConstraint) = 1
minimal_npairs(::HexagonalConstraint) = 2
minimal_npairs(::TrigonalConstraint) = 2
minimal_npairs(::TetragonalConstraint) = 2
minimal_npairs(::OrthorhombicConstraint) = 3
minimal_npairs(::MonoclinicConstraint) = 5
minimal_npairs(::TriclinicConstraint) = 6
