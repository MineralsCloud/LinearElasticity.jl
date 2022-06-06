using Crystallography:
    CrystalSystem,
    Cubic,
    Hexagonal,
    Trigonal,
    Tetragonal,
    Orthorhombic,
    Monoclinic,
    Triclinic

export solve_elastic_constants, solve_stiffnesses, solve_compliances

function recombine_strains(::Cubic, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    return [  # 6×3 matrix
        ϵ₁ ϵ₂+ϵ₃ 0
        ϵ₂ ϵ₁+ϵ₃ 0
        ϵ₃ ϵ₁+ϵ₂ 0
        0 0 ϵ₄
        0 0 ϵ₅
        0 0 ϵ₆
    ]
end
function recombine_strains(::Tetragonal, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    # Tetragonal (I) class (c₁₆ = 0) is a special case of tetragonal (II) class
    return [  # 6×7 matrix
        ϵ₁ 0 ϵ₂ ϵ₃ ϵ₆ 0 0
        ϵ₂ 0 ϵ₁ ϵ₃ -ϵ₆ 0 0
        0 ϵ₃ 0 ϵ₁+ϵ₂ 0 0 0
        0 0 0 0 0 ϵ₄ 0
        0 0 0 0 0 ϵ₅ 0
        0 0 0 0 ϵ₁-ϵ₂ 0 ϵ₆
    ]
end
function recombine_strains(::Orthorhombic, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    return [  # 6×9 matrix
        ϵ₁ 0 0 ϵ₂ ϵ₃ 0 0 0 0
        0 ϵ₂ 0 ϵ₁ 0 ϵ₃ 0 0 0
        0 0 ϵ₃ 0 ϵ₁ ϵ₂ 0 0 0
        0 0 0 0 0 0 ϵ₄ 0 0
        0 0 0 0 0 0 0 ϵ₅ 0
        0 0 0 0 0 0 0 0 ϵ₆
    ]
end
function recombine_strains(::Hexagonal, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    return [  # 6×5 matrix
        ϵ₁ 0 ϵ₂ ϵ₃ 0
        ϵ₂ 0 ϵ₁ ϵ₃ 0
        0 ϵ₃ 0 ϵ₁+ϵ₂ 0
        0 0 0 0 ϵ₄
        0 0 0 0 ϵ₅
        ϵ₆/2 0 -ϵ₆/2 0 0
    ]
end
function recombine_strains(::Trigonal, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    # Rhombohedral (I) class (c₁₅ = 0) is a special case of rhombohedral (II) class
    return [  # 6×7 matrix
        ϵ₁ 0 ϵ₂ ϵ₃ 0 ϵ₄ ϵ₅
        ϵ₂ 0 ϵ₁ ϵ₃ 0 -ϵ₄ -ϵ₅
        0 ϵ₃ 0 ϵ₁+ϵ₂ 0 0 0
        0 0 0 0 ϵ₄ ϵ₁-ϵ₂ -ϵ₆
        0 0 0 0 ϵ₅ ϵ₆ ϵ₁-ϵ₂
        ϵ₆/2 0 -ϵ₆/2 0 0 ϵ₅ -ϵ₄
    ]
end
function recombine_strains(::Monoclinic, strain::EngineeringStrain)  # Only standard orientation (diad ∥ x₂) is implemented
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    return [  # 6×13 matrix
        ϵ₁ 0 0 ϵ₂ ϵ₃ 0 0 0 0 ϵ₅ 0 0 0
        0 ϵ₂ 0 ϵ₁ 0 ϵ₃ 0 0 0 0 ϵ₅ 0 0
        0 0 ϵ₃ 0 ϵ₁ ϵ₂ 0 0 0 0 0 ϵ₅ 0
        0 0 0 0 0 0 ϵ₄ 0 0 0 0 0 ϵ₆
        0 0 0 0 0 0 0 ϵ₅ 0 ϵ₁ ϵ₂ ϵ₃ 0
        0 0 0 0 0 0 0 0 ϵ₆ 0 0 0 ϵ₄
    ]
end
function recombine_strains(::Triclinic, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    return [  # 6×21 matrix
        ϵ₁ ϵ₂ ϵ₃ ϵ₄ ϵ₅ ϵ₆ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 ϵ₁ 0 0 0 0 ϵ₂ ϵ₃ ϵ₄ ϵ₅ ϵ₆ 0 0 0 0 0 0 0 0 0 0
        0 0 ϵ₁ 0 0 0 0 ϵ₂ 0 0 0 ϵ₃ ϵ₄ ϵ₅ ϵ₆ 0 0 0 0 0 0
        0 0 0 ϵ₁ 0 0 0 0 ϵ₂ 0 0 0 ϵ₃ 0 0 ϵ₄ ϵ₅ ϵ₆ 0 0 0
        0 0 0 0 ϵ₁ 0 0 0 0 ϵ₂ 0 0 0 ϵ₃ 0 0 ϵ₄ 0 ϵ₅ ϵ₆ 0
        0 0 0 0 0 ϵ₁ 0 0 0 0 ϵ₂ 0 0 0 ϵ₃ 0 0 ϵ₄ 0 ϵ₅ ϵ₆
    ]
end
recombine_strains(system::CrystalSystem, strains::AbstractVector{<:EngineeringStrain}) =
    vcat((recombine_strains(system, strain) for strain in strains)...)

function reconstruct_cᵢⱼ(::Cubic, coefficients)
    c₁₁, c₁₂, c₄₄ = coefficients
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        c₁₁,
        c₁₂,
        c₁₂,
        𝟎,
        𝟎,
        𝟎,
        c₁₁,
        c₁₂,
        𝟎,
        𝟎,
        𝟎,
        c₁₁,
        𝟎,
        𝟎,
        𝟎,
        c₄₄,
        𝟎,
        𝟎,
        c₄₄,
        𝟎,
        c₄₄,
    )
end
function reconstruct_cᵢⱼ(::Tetragonal, coefficients)
    c₁₁, c₃₃, c₁₂, c₁₃, c₁₆, c₄₄, c₆₆ = coefficients
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        c₁₁,
        c₁₂,
        c₁₃,
        𝟎,
        𝟎,
        c₁₆,
        c₁₁,
        c₁₃,
        𝟎,
        𝟎,
        -c₁₆,
        c₃₃,
        𝟎,
        𝟎,
        𝟎,
        c₄₄,
        𝟎,
        𝟎,
        c₄₄,
        𝟎,
        c₆₆,
    )
end
function reconstruct_cᵢⱼ(::Orthorhombic, coefficients)
    c₁₁, c₂₂, c₃₃, c₁₂, c₁₃, c₂₃, c₄₄, c₅₅, c₆₆ = coefficients
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        c₁₁,
        c₁₂,
        c₁₃,
        𝟎,
        𝟎,
        𝟎,
        c₂₂,
        c₂₃,
        𝟎,
        𝟎,
        𝟎,
        c₃₃,
        𝟎,
        𝟎,
        𝟎,
        c₄₄,
        𝟎,
        𝟎,
        c₅₅,
        𝟎,
        c₆₆,
    )
end
function reconstruct_cᵢⱼ(::Hexagonal, coefficients)
    c₁₁, c₃₃, c₁₂, c₁₃, c₄₄ = coefficients
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        c₁₁,
        c₁₂,
        c₁₃,
        𝟎,
        𝟎,
        𝟎,
        c₁₁,
        c₁₃,
        𝟎,
        𝟎,
        𝟎,
        c₃₃,
        𝟎,
        𝟎,
        𝟎,
        c₄₄,
        𝟎,
        𝟎,
        c₄₄,
        𝟎,
        (c₁₁ - c₁₂) / 2,
    )
end
function reconstruct_cᵢⱼ(::Trigonal, coefficients)
    c₁₁, c₃₃, c₁₂, c₁₃, c₄₄, c₁₄, c₁₅ = coefficients
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        c₁₁,
        c₁₂,
        c₁₃,
        c₁₄,
        c₁₅,
        𝟎,
        c₁₁,
        c₁₃,
        -c₁₄,
        -c₁₅,
        𝟎,
        c₃₃,
        𝟎,
        𝟎,
        𝟎,
        c₄₄,
        𝟎,
        -c₁₅,
        c₄₄,
        c₁₄,
        (c₁₁ - c₁₂) / 2,
    )
end
function reconstruct_cᵢⱼ(::Monoclinic, coefficients)
    c₁₁, c₂₂, c₃₃, c₁₂, c₁₃, c₂₃, c₄₄, c₅₅, c₆₆, c₁₅, c₂₅, c₃₅, c₄₆ = coefficients
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        c₁₁,
        c₁₂,
        c₁₃,
        𝟎,
        c₁₅,
        𝟎,
        c₂₂,
        c₂₃,
        𝟎,
        c₂₅,
        𝟎,
        c₃₃,
        𝟎,
        c₃₅,
        𝟎,
        c₄₄,
        𝟎,
        c₄₆,
        c₅₅,
        𝟎,
        c₆₆,
    )
end
reconstruct_cᵢⱼ(::Triclinic, coefficients) = StiffnessMatrix(coefficients...)

function recombine_stresses(::Cubic, stress::EngineeringStress)
    σ₁, σ₂, σ₃, σ₄, σ₅, σ₆ = stress
    return [  # 6×3 matrix
        σ₁ σ₂+σ₃ 0
        σ₂ σ₁+σ₃ 0
        σ₃ σ₁+σ₂ 0
        0 0 σ₄
        0 0 σ₅
        0 0 σ₆
    ]
end
recombine_stresses(system::CrystalSystem, stresses::AbstractVector{<:EngineeringStress}) =
    vcat((recombine_stresses(system, stress) for stress in stresses)...)

function reconstruct_sᵢⱼ(::Cubic, coefficients)
    s₁₁, s₁₂, s₄₄ = coefficients
    𝟎 = zero(s₁₁)
    return StiffnessMatrix(
        s₁₁,
        s₁₂,
        s₁₂,
        𝟎,
        𝟎,
        𝟎,
        s₁₁,
        s₁₂,
        𝟎,
        𝟎,
        𝟎,
        s₁₁,
        𝟎,
        𝟎,
        𝟎,
        s₄₄,
        𝟎,
        𝟎,
        s₄₄,
        𝟎,
        s₄₄,
    )
end

function solve_elastic_constants(
    system::CrystalSystem,
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    if length(strains) != length(stresses)
        throw(
            DimensionMismatch(
                "the number of strains and the number of stresses must match!",
            ),
        )
    end
    n = minimal_npairs(system)
    if length(strains) < n
        throw(ArgumentError("the number of strains/stresses must be at least $n."))
    end
    σ = vcat(stresses...)  # Length 6n vector, n = length(strains) = length(stresses)
    ε = recombine_strains(system, strains)  # Size 6n×N matrix, N = # independent coefficients
    coefficients = ε \ σ  # Length N vector
    return reconstruct_cᵢⱼ(system, coefficients)
end
function solve_elastic_constants(
    system::CrystalSystem,
    stresses::AbstractVector{<:EngineeringStress},
    strains::AbstractVector{<:EngineeringStrain},
)
    if length(strains) != length(stresses)
        throw(
            DimensionMismatch(
                "the number of strains and the number of stresses must match!",
            ),
        )
    end
    n = minimal_npairs(system)
    if length(strains) < n
        throw(ArgumentError("the number of strains/stresses must be at least $n."))
    end
    ε = vcat(strains...)
    σ = recombine_stresses(system, stresses)
    coefficients = σ \ ε
    return reconstruct_sᵢⱼ(system, coefficients)
end
function solve_elastic_constants(
    system::CrystalSystem,
    strains::AbstractVector{<:TensorStrain},
    stresses::AbstractVector{<:TensorStress},
)
    cᵢⱼ = solve_elastic_constants(
        system,
        EngineeringStrain.(strains),
        EngineeringStress.(stresses),
    )
    return StiffnessTensor(cᵢⱼ)
end
function solve_elastic_constants(
    system::CrystalSystem,
    stresses::AbstractVector{<:TensorStress},
    strains::AbstractVector{<:TensorStrain},
)
    sᵢⱼ = solve_elastic_constants(
        system,
        EngineeringStrain.(strains),
        EngineeringStress.(stresses),
    )
    return ComplianceTensor(sᵢⱼ)
end
solve_elastic_constants(strains_or_stresses, stresses_or_strains) =
    solve_elastic_constants(Triclinic(), strains_or_stresses, stresses_or_strains)

solve_stiffnesses(
    system::CrystalSystem,
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
) = solve_elastic_constants(system, strains, stresses)
solve_stiffnesses(
    system::CrystalSystem,
    strains::AbstractVector{<:TensorStrain},
    stresses::AbstractVector{<:TensorStress},
) = solve_elastic_constants(system, strains, stresses)
solve_stiffnesses(strains, stresses) = solve_stiffnesses(Triclinic(), strains, stresses)

solve_compliances(
    system::CrystalSystem,
    stresses::AbstractVector{<:EngineeringStress},
    strains::AbstractVector{<:EngineeringStrain},
) = solve_elastic_constants(system, stresses, strains)
solve_compliances(
    system::CrystalSystem,
    stresses::AbstractVector{<:TensorStress},
    strains::AbstractVector{<:TensorStrain},
) = solve_elastic_constants(system, stresses, strains)
solve_compliances(stresses, strains) = solve_compliances(Triclinic(), stresses, strains)

minimal_npairs(::Cubic) = 1
minimal_npairs(::Hexagonal) = 2
minimal_npairs(::Trigonal) = 2
minimal_npairs(::Tetragonal) = 2
minimal_npairs(::Orthorhombic) = 3
minimal_npairs(::Monoclinic) = 5
minimal_npairs(::Triclinic) = 6
