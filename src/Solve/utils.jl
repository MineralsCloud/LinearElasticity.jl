function make_linear_operator(ϵ::EngineeringStrain, ::Cubic)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = ϵ
    return [  # 6×3 matrix
        ϵ₁ ϵ₂+ϵ₃ 0
        ϵ₂ ϵ₁+ϵ₃ 0
        ϵ₃ ϵ₁+ϵ₂ 0
        0 0 ϵ₄
        0 0 ϵ₅
        0 0 ϵ₆
    ]
end
function make_linear_operator(ϵ::EngineeringStrain, ::Tetragonal)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = ϵ
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
function make_linear_operator(ϵ::EngineeringStrain, ::Orthorhombic)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = ϵ
    return [  # 6×9 matrix
        ϵ₁ 0 0 ϵ₂ ϵ₃ 0 0 0 0
        0 ϵ₂ 0 ϵ₁ 0 ϵ₃ 0 0 0
        0 0 ϵ₃ 0 ϵ₁ ϵ₂ 0 0 0
        0 0 0 0 0 0 ϵ₄ 0 0
        0 0 0 0 0 0 0 ϵ₅ 0
        0 0 0 0 0 0 0 0 ϵ₆
    ]
end
function make_linear_operator(ϵ::EngineeringStrain, ::Hexagonal)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = ϵ
    return [  # 6×5 matrix
        ϵ₁ 0 ϵ₂ ϵ₃ 0
        ϵ₂ 0 ϵ₁ ϵ₃ 0
        0 ϵ₃ 0 ϵ₁+ϵ₂ 0
        0 0 0 0 ϵ₄
        0 0 0 0 ϵ₅
        ϵ₆/2 0 -ϵ₆/2 0 0
    ]
end
function make_linear_operator(ϵ::EngineeringStrain, ::Trigonal)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = ϵ
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
function make_linear_operator(ϵ::EngineeringStrain, ::Monoclinic)  # Only standard orientation (diad ∥ x₂) is implemented
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = ϵ
    return [  # 6×13 matrix
        ϵ₁ 0 0 ϵ₂ ϵ₃ 0 0 0 0 ϵ₅ 0 0 0
        0 ϵ₂ 0 ϵ₁ 0 ϵ₃ 0 0 0 0 ϵ₅ 0 0
        0 0 ϵ₃ 0 ϵ₁ ϵ₂ 0 0 0 0 0 ϵ₅ 0
        0 0 0 0 0 0 ϵ₄ 0 0 0 0 0 ϵ₆
        0 0 0 0 0 0 0 ϵ₅ 0 ϵ₁ ϵ₂ ϵ₃ 0
        0 0 0 0 0 0 0 0 ϵ₆ 0 0 0 ϵ₄
    ]
end
function make_linear_operator(ϵ::EngineeringStrain, ::Triclinic)
    ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆ = ϵ
    return [  # 6×21 matrix
        ϵ₁ ϵ₂ ϵ₃ ϵ₄ ϵ₅ ϵ₆ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 ϵ₁ 0 0 0 0 ϵ₂ ϵ₃ ϵ₄ ϵ₅ ϵ₆ 0 0 0 0 0 0 0 0 0 0
        0 0 ϵ₁ 0 0 0 0 ϵ₂ 0 0 0 ϵ₃ ϵ₄ ϵ₅ ϵ₆ 0 0 0 0 0 0
        0 0 0 ϵ₁ 0 0 0 0 ϵ₂ 0 0 0 ϵ₃ 0 0 ϵ₄ ϵ₅ ϵ₆ 0 0 0
        0 0 0 0 ϵ₁ 0 0 0 0 ϵ₂ 0 0 0 ϵ₃ 0 0 ϵ₄ 0 ϵ₅ ϵ₆ 0
        0 0 0 0 0 ϵ₁ 0 0 0 0 ϵ₂ 0 0 0 ϵ₃ 0 0 ϵ₄ 0 ϵ₅ ϵ₆
    ]
end
make_linear_operator(
    𝛜::AbstractVector{<:EngineeringStrain}, constraint::SymmetryConstraint
) = vcat((make_linear_operator(ϵ, constraint) for ϵ in 𝛜)...)
function make_linear_operator(σ::EngineeringStress, ::Cubic)
    σ₁, σ₂, σ₃, σ₄, σ₅, σ₆ = σ
    return [  # 6×3 matrix
        σ₁ σ₂+σ₃ 0
        σ₂ σ₁+σ₃ 0
        σ₃ σ₁+σ₂ 0
        0 0 σ₄
        0 0 σ₅
        0 0 σ₆
    ]
end
make_linear_operator(
    𝛔::AbstractVector{<:EngineeringStress}, constraint::SymmetryConstraint
) = vcat((make_linear_operator(σ, constraint) for σ in 𝛔)...)

function construct_cᵢⱼ(𝐜, ::Cubic)
    𝟎, c₁₁, c₁₂, c₄₄ = _promote_with_zero(𝐜)
    return StiffnessMatrix(
        c₁₁, c₁₂, c₁₂, 𝟎, 𝟎, 𝟎, c₁₁, c₁₂, 𝟎, 𝟎, 𝟎, c₁₁, 𝟎, 𝟎, 𝟎, c₄₄, 𝟎, 𝟎, c₄₄, 𝟎, c₄₄
    )
end
function construct_cᵢⱼ(𝐜, ::Tetragonal)
    𝟎, c₁₁, c₃₃, c₁₂, c₁₃, c₁₆, c₄₄, c₆₆ = _promote_with_zero(𝐜)
    return StiffnessMatrix(
        c₁₁, c₁₂, c₁₃, 𝟎, 𝟎, c₁₆, c₁₁, c₁₃, 𝟎, 𝟎, -c₁₆, c₃₃, 𝟎, 𝟎, 𝟎, c₄₄, 𝟎, 𝟎, c₄₄, 𝟎, c₆₆
    )
end
function construct_cᵢⱼ(𝐜, ::Orthorhombic)
    𝟎, c₁₁, c₂₂, c₃₃, c₁₂, c₁₃, c₂₃, c₄₄, c₅₅, c₆₆ = _promote_with_zero(𝐜)
    return StiffnessMatrix(
        c₁₁, c₁₂, c₁₃, 𝟎, 𝟎, 𝟎, c₂₂, c₂₃, 𝟎, 𝟎, 𝟎, c₃₃, 𝟎, 𝟎, 𝟎, c₄₄, 𝟎, 𝟎, c₅₅, 𝟎, c₆₆
    )
end
function construct_cᵢⱼ(𝐜, ::Hexagonal)
    𝟎, c₁₁, c₃₃, c₁₂, c₁₃, c₄₄ = _promote_with_zero(𝐜)
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
function construct_cᵢⱼ(𝐜, ::Trigonal)
    𝟎, c₁₁, c₃₃, c₁₂, c₁₃, c₄₄, c₁₄, c₁₅ = _promote_with_zero(𝐜)
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
function construct_cᵢⱼ(𝐜, ::Monoclinic)
    𝟎, c₁₁, c₂₂, c₃₃, c₁₂, c₁₃, c₂₃, c₄₄, c₅₅, c₆₆, c₁₅, c₂₅, c₃₅, c₄₆ = _promote_with_zero(
        𝐜
    )
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
construct_cᵢⱼ(𝐜, ::Triclinic) = StiffnessMatrix(𝐜...)

function construct_sᵢⱼ(𝐬, ::Cubic)
    𝟎, s₁₁, s₁₂, s₄₄ = _promote_with_zero(𝐬)
    return ComplianceMatrix(
        s₁₁, s₁₂, s₁₂, 𝟎, 𝟎, 𝟎, s₁₁, s₁₂, 𝟎, 𝟎, 𝟎, s₁₁, 𝟎, 𝟎, 𝟎, s₄₄, 𝟎, 𝟎, s₄₄, 𝟎, s₄₄
    )
end

function _promote_with_zero(xs)
    T = Base.promote_typeof(xs...)
    𝟎 = zero(T)
    return 𝟎, (convert(T, x) for x in xs)...
end
