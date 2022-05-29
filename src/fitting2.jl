using Crystallography:
    CrystalSystem,
    Cubic,
    Hexagonal,
    Trigonal,
    Tetragonal,
    Orthorhombic,
    Monoclinic,
    Triclinic
using LinearAlgebra: Symmetric, dot

export fit_elastic_constant

function form_matrix(::Cubic, strain::EngineeringStrain)
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
function form_matrix(::Tetragonal, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    return [  # 6×6 matrix
        ϵ₁ 0 ϵ₂ ϵ₃ 0 0
        ϵ₂ 0 ϵ₁ ϵ₃ 0 0
        0 ϵ₃ 0 ϵ₁+ϵ₂ 0 0
        0 0 0 0 0 ϵ₆
        0 0 0 0 ϵ₅ 0
        0 0 0 0 ϵ₄ 0
    ]
end
function form_matrix(::Orthorhombic, strain::EngineeringStrain)
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
function form_matrix(::Hexagonal, strain::EngineeringStrain)
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
function form_matrix(::Trigonal, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    return [  # 6×6 matrix
        ϵ₁ 0 ϵ₂ ϵ₃ 0 ϵ₅
        ϵ₂ 0 ϵ₁ ϵ₃ 0 -ϵ₅
        0 ϵ₃ 0 ϵ₁+ϵ₂ 0 0
        0 0 0 0 ϵ₄ -2ϵ₆
        0 0 0 0 ϵ₅ 2(ϵ₁-ϵ₂)
        ϵ₆ 0 -ϵ₆ 0 0 -2ϵ₄
    ]
end
function form_matrix(::Monoclinic, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    γ = ϵ₆ / 2
    return [  # 6×13 matrix
        ϵ₁ 0 0 ϵ₂ ϵ₃ 0 0 0 0 γ 0 0 0
        0 ϵ₂ 0 ϵ₁ 0 ϵ₃ 0 0 0 0 γ 0 0
        0 0 ϵ₃ 0 ϵ₁ ϵ₂ 0 0 0 0 0 γ 0
        0 0 0 0 0 0 ϵ₄ 0 0 0 0 0 ϵ₅/2
        0 0 0 0 0 0 0 ϵ₅ 0 0 0 0 ϵ₄/2
        0 0 0 0 0 0 0 0 ϵ₆ ϵ₁ ϵ₂ ϵ₃ 0
    ]
end
function form_matrix(::Triclinic, strain::EngineeringStrain)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = strain
    α, β, γ = ϵ₄ / 2, ϵ₅ / 2, ϵ₆ / 2
    return [  # 6×21 matrix
        ϵ₁ 0 0 ϵ₂ ϵ₃ 0 0 0 0 γ 0 0 0 0 α β 0 0
        0 ϵ₂ 0 ϵ₁ 0 ϵ₃ 0 0 0 0 γ 0 0 0 0 0 β 0
        0 0 ϵ₃ 0 ϵ₁ ϵ₂ 0 0 0 0 0 γ 0 0 0 0 0 0
        0 0 0 0 0 0 ϵ₄ 0 0 0 0 0 γ 0 ϵ₁ 0 0 β
        0 0 0 0 0 0 0 ϵ₅ 0 0 0 0 0 γ 0 ϵ₁ ϵ₂ α
        0 0 0 0 0 0 0 0 ϵ₆ ϵ₁ ϵ₂ ϵ₃ α β 0 0 0 0
    ]
end
form_matrix(system::CrystalSystem, strains::AbstractVector{<:EngineeringStrain}) =
    vcat(form_matrix(system, strain) for strain in strains)

function reorder_cᵢⱼ(::Cubic, cᵢⱼ)
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
function reorder_cᵢⱼ(::Tetragonal, cᵢⱼ)
    c₁₁, c₃₃, c₁₂, c₁₃, c₄₄, c₁₄ = cᵢⱼ
    𝟎 = zero(c₁₁)
    return StiffnessMatrix([])
end
function reorder_cᵢⱼ(::Orthorhombic, cᵢⱼ)
    c₁₁, c₂₂, c₃₃, c₁₂, c₁₃, c₂₃, c₄₄, c₅₅, c₆₆ = cᵢⱼ
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        [
            c₁₁ c₁₂ c₁₃ 𝟎 𝟎 𝟎
            c₁₂ c₂₂ c₂₃ 𝟎 𝟎 𝟎
            c₁₃ c₂₃ c₃₃ 𝟎 𝟎 𝟎
            𝟎 𝟎 𝟎 c₄₄ 𝟎 𝟎
            𝟎 𝟎 𝟎 𝟎 c₅₅ 𝟎
            𝟎 𝟎 𝟎 𝟎 𝟎 c₆₆
        ],
        c₁₁ c₁₂ c₁₃ 𝟎 𝟎 𝟎 c₂₂ c₂₃ 𝟎 𝟎 𝟎 c₃₃ 𝟎 𝟎 𝟎 c₄₄ 𝟎 𝟎 c₅₅ 𝟎
    )
end
function reorder_cᵢⱼ(::Hexagonal, cᵢⱼ)
    c₁₁, c₃₃, c₁₂, c₁₃, c₄₄ = cᵢⱼ
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        [
            c₁₁ c₁₂ c₁₃ 𝟎 𝟎 𝟎
            c₁₂ c₁₁ c₁₃ 𝟎 𝟎 𝟎
            c₁₃ c₁₃ c₃₃ 𝟎 𝟎 𝟎
            𝟎 𝟎 𝟎 c₄₄ 𝟎 𝟎
            𝟎 𝟎 𝟎 𝟎 c₄₄ 𝟎
            𝟎 𝟎 𝟎 𝟎 𝟎 (c₁₁-c₁₂)/2
        ],
    )
end
function reorder_cᵢⱼ(::Trigonal, cᵢⱼ)
    c₁₁, c₃₃, c₁₂, c₁₃, c₄₄, c₁₄ = cᵢⱼ
    𝟎 = zero(c₁₁)
    return StiffnessMatrix([])
end
function reorder_cᵢⱼ(::Monoclinic, cᵢⱼ)
    c₁₁, c₂₂, c₃₃, c₁₂, c₁₃, c₂₃, c₄₄, c₅₅, c₆₆, c₁₆, c₂₆, c₃₆, c₄₅ = cᵢⱼ
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        [
            c₁₁ c₁₂ c₁₃ 𝟎 𝟎 c₁₆
            c₁₂ c₂₂ c₂₃ 𝟎 𝟎 c₂₆
            c₁₃ c₂₃ c₃₃ 𝟎 𝟎 c₃₆
            𝟎 𝟎 𝟎 c₄₄ c₄₅ 𝟎
            𝟎 𝟎 𝟎 𝟎 c₅₅ 𝟎
            𝟎 𝟎 𝟎 𝟎 𝟎 c₆₆
        ],
    )
end
function reorder_cᵢⱼ(::Triclinic, cᵢⱼ)
    c₁₁,
    c₂₂,
    c₃₃,
    c₁₂,
    c₁₃,
    c₂₃,
    c₄₄,
    c₅₅,
    c₆₆,
    c₁₆,
    c₂₆,
    c₃₆,
    c₄₆,
    c₅₆,
    c₁₄,
    c₁₅,
    c₂₅,
    c₄₅ = cᵢⱼ
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        [
            c₁₁ c₁₂ c₁₃ c₁₄ 𝟎 c₁₆
            c₁₂ c₂₂ c₂₃ 𝟎 𝟎 c₂₆
            c₁₃ c₂₃ c₃₃ 𝟎 𝟎 c₃₆
            𝟎 𝟎 𝟎 c₄₄ c₄₅ 𝟎
            𝟎 𝟎 𝟎 𝟎 c₅₅ 𝟎
            𝟎 𝟎 𝟎 𝟎 𝟎 c₆₆
        ],
    )
end

function fit_elastic_constant(
    system::CrystalSystem,
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    σ = vcat(stresses...)  # Length 6n vector, n = length(strains) = length(stresses)
    ε = form_matrix(system, strains)  # Size 6n×N matrix, N = # independent coefficients
    cᵢⱼ = ε \ σ  # Length N vector
    return reorder_cᵢⱼ(system, cᵢⱼ)
end
