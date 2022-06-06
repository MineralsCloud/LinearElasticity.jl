using Crystallography:
    CrystalSystem,
    Cubic,
    Hexagonal,
    Trigonal,
    Tetragonal,
    Orthorhombic,
    Monoclinic,
    Triclinic

export solve_elastic_matrix

function construct_strain_matrix(::Cubic, strain::EngineeringStrain)
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
function construct_strain_matrix(::Tetragonal, strain::EngineeringStrain)
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
function construct_strain_matrix(::Orthorhombic, strain::EngineeringStrain)
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
function construct_strain_matrix(::Hexagonal, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    return [  # 6×5 matrix
        ϵ₁ 0 ϵ₂ ϵ₃ 0
        ϵ₂ 0 ϵ₁ ϵ₃ 0
        0 ϵ₃ 0 ϵ₁+ϵ₂ 0
        0 0 0 0 ϵ₄
        0 0 0 0 ϵ₅
        ϵ₆ 0 -ϵ₆ 0 0
    ]
end
function construct_strain_matrix(::Trigonal, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    # Rhombohedral (I) class (c₁₅ = 0) is a special case of rhombohedral (II) class
    return [  # 6×6 matrix
        ϵ₁ 0 ϵ₂ ϵ₃ 0 ϵ₄ ϵ₅
        ϵ₂ 0 ϵ₁ ϵ₃ 0 -ϵ₄ -ϵ₅
        0 ϵ₃ 0 ϵ₁+ϵ₂ 0 0 0
        0 0 0 0 ϵ₄ ϵ₁-ϵ₂ -ϵ₆
        0 0 0 0 ϵ₅ ϵ₆ ϵ₁-ϵ₂
        ϵ₆/2 0 -ϵ₆/2 0 0 ϵ₅ -ϵ₄
    ]
end
function construct_strain_matrix(::Monoclinic, strain::EngineeringStrain)  # Only standard orientation (diad ∥ x₂) is implemented
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
function construct_strain_matrix(::Triclinic, strain::EngineeringStrain)
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
construct_strain_matrix(
    system::CrystalSystem,
    strains::AbstractVector{<:EngineeringStrain},
) = vcat((construct_strain_matrix(system, strain) for strain in strains)...)

function reconstruct_cᵢⱼ(::Cubic, cᵢⱼ)
    c₁₁, c₁₂, c₄₄ = cᵢⱼ
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
function reconstruct_cᵢⱼ(::Tetragonal, cᵢⱼ)
    c₁₁, c₃₃, c₁₂, c₁₃, c₁₆, c₄₄, c₆₆ = cᵢⱼ
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
function reconstruct_cᵢⱼ(::Orthorhombic, cᵢⱼ)
    c₁₁, c₂₂, c₃₃, c₁₂, c₁₃, c₂₃, c₄₄, c₅₅, c₆₆ = cᵢⱼ
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
function reconstruct_cᵢⱼ(::Hexagonal, cᵢⱼ)
    c₁₁, c₃₃, c₁₂, c₁₃, c₄₄ = cᵢⱼ
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
function reconstruct_cᵢⱼ(::Trigonal, cᵢⱼ)
    c₁₁, c₃₃, c₁₂, c₁₃, c₄₄, c₁₄, c₁₅ = cᵢⱼ
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
function reconstruct_cᵢⱼ(::Monoclinic, cᵢⱼ)
    c₁₁, c₂₂, c₃₃, c₁₂, c₁₃, c₂₃, c₄₄, c₅₅, c₆₆, c₁₅, c₂₅, c₃₅, c₄₆ = cᵢⱼ
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
reconstruct_cᵢⱼ(::Triclinic, cᵢⱼ) = StiffnessMatrix(cᵢⱼ...)

function construct_stress_matrix(::Cubic, stress::EngineeringStress)
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
construct_stress_matrix(
    system::CrystalSystem,
    stresses::AbstractVector{<:EngineeringStrain},
) = vcat(construct_stress_matrix(system, stress) for stress in stresses)

function reconstruct_sᵢⱼ(::Cubic, sᵢⱼ)
    s₁₁, s₁₂, s₄₄ = sᵢⱼ
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

function solve_elastic_matrix(
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
    n = minimal_ulics(system)
    if length(strains) < n
        throw(ArgumentError("the number of strains/stresses must be at least $n."))
    end
    σ = vcat(stresses...)  # Length 6n vector, n = length(strains) = length(stresses)
    ε = construct_strain_matrix(system, strains)  # Size 6n×N matrix, N = # independent coefficients
    cᵢⱼ = ε \ σ  # Length N vector
    return reconstruct_cᵢⱼ(system, cᵢⱼ)
end
function solve_elastic_matrix(
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
    n = minimal_ulics(system)
    if length(strains) < n
        throw(ArgumentError("the number of strains/stresses must be at least $n."))
    end
    ε = vcat(strains...)
    σ = construct_stress_matrix(system, stresses)
    sᵢⱼ = σ \ ε
    return reconstruct_sᵢⱼ(system, sᵢⱼ)
end
solve_elastic_matrix(
    system::CrystalSystem,
    strains::AbstractVector{<:TensorStrain},
    stresses::AbstractVector{<:TensorStress},
) = solve_elastic_matrix(system, EngineeringStrain.(strains), EngineeringStress.(stresses))
solve_elastic_matrix(
    system::CrystalSystem,
    stresses::AbstractVector{<:TensorStress},
    strains::AbstractVector{<:TensorStrain},
) = solve_elastic_matrix(system, EngineeringStress.(stresses), EngineeringStrain.(strains))
solve_elastic_matrix(strains_or_stresses, stresses_or_strains) =
    solve_elastic_matrix(Triclinic(), strains_or_stresses, stresses_or_strains)

minimal_ulics(::Cubic) = 1
minimal_ulics(::Hexagonal) = 2
minimal_ulics(::Trigonal) = 2
minimal_ulics(::Tetragonal) = 2
minimal_ulics(::Orthorhombic) = 3
minimal_ulics(::Monoclinic) = 5
minimal_ulics(::Triclinic) = 6
