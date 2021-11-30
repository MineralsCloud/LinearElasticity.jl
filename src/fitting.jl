using CrystallographyBase: CrystalSystem, Cubic, Hexagonal, Orthorhombic
using LinearAlgebra: Symmetric, dot

export ElasticConstantFitter

struct ElasticConstantFitter{T<:CrystalSystem}
    system::T
end

function (::ElasticConstantFitter{Orthorhombic})(
    ϵs::AbstractVector{<:EngineeringStrain},
    σs::AbstractVector{<:EngineeringStress},
)
    @assert length(ϵs) == length(σs) >= 3
    ϵ, σ = hcat(ϵs...), hcat(σs...)
    c = transpose(ϵ[1:3, :]) \ transpose(σ[1:3, :]) |> Symmetric
    Bᵀ = transpose(ϵ[4:6, :])
    c₄₄ = (Bᵀ\σ[4, :])[1]
    c₅₅ = (Bᵀ\σ[5, :])[2]
    c₆₆ = (Bᵀ\σ[6, :])[3]
    𝟎 = zero(c[1, 1])
    mat = vcat(
        hcat(c, fill(𝟎, 3, 3)),
        [
            𝟎 𝟎 𝟎 c₄₄ 𝟎 𝟎
            𝟎 𝟎 𝟎 𝟎 c₅₅ 𝟎
            𝟎 𝟎 𝟎 𝟎 𝟎 c₆₆
        ],
    )
    return StiffnessMatrix(mat)
end
function (::ElasticConstantFitter{Hexagonal})(
    ϵs::AbstractVector{<:EngineeringStrain},
    σs::AbstractVector{<:EngineeringStress},
)
    @assert length(ϵs) == length(σs) >= 2
    ϵ₁, ϵ₂, ϵ₃ = ϵs[1][1:3]
    ϵ₁′, ϵ₂′, ϵ₃′ = ϵs[2][1:3]
    Aᵀ = [
        ϵ₁ ϵ₂ 0 ϵ₁′ ϵ₂′ 0
        ϵ₂ ϵ₁ 0 ϵ₂′ ϵ₁′ 0
        ϵ₃ ϵ₃ ϵ₁+ϵ₂ ϵ₃′ ϵ₃′ ϵ₁′+ϵ₂′
        0 0 ϵ₃ 0 0 ϵ₃′
    ]
    c₁₁, c₁₂, c₁₃, c₃₃ = inv(Aᵀ * transpose(Aᵀ)) * Aᵀ * append!(σs[1][1:3], σs[2][1:3])
    c₄₄ = σs[1][4] / ϵs[1][4]
    c₆₆ = σs[1][6] / ϵs[1][6]
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        [
            c₁₁ c₁₂ c₁₃ 𝟎 𝟎 𝟎
            c₁₂ c₁₁ c₁₃ 𝟎 𝟎 𝟎
            c₁₃ c₁₃ c₃₃ 𝟎 𝟎 𝟎
            𝟎 𝟎 𝟎 c₄₄ 𝟎 𝟎
            𝟎 𝟎 𝟎 𝟎 c₄₄ 𝟎
            𝟎 𝟎 𝟎 𝟎 𝟎 c₆₆
        ],
    )
end
function (::ElasticConstantFitter{Cubic})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    ϵ, σ = first(strains), first(stresses)
    ϵ₁, ϵ₂, ϵ₃ = ϵ[1:3]
    Aᵀ = [
        ϵ₁ ϵ₂ ϵ₃
        ϵ₂+ϵ₃ ϵ₁+ϵ₃ ϵ₂+ϵ₁
    ]
    c₁₁, c₁₂ = inv(Aᵀ * transpose(Aᵀ)) * Aᵀ * σ[1:3]  # If 𝐴 is well-conditioned, using the normal equations is around as accurate as other methods and is also the fastest. https://math.stackexchange.com/a/3252377/115512
    c₄₄ = dot(ϵ[4:6], σ[4:6]) / sum(abs2, ϵ[4:6])  # B = ϵ[4:6], c₄₄ = inv(Bᵀ * B) * Bᵀ * σ[4:6]
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        [
            c₁₁ c₁₂ c₁₂ 𝟎 𝟎 𝟎
            c₁₂ c₁₁ c₁₂ 𝟎 𝟎 𝟎
            c₁₂ c₁₂ c₁₁ 𝟎 𝟎 𝟎
            𝟎 𝟎 𝟎 c₄₄ 𝟎 𝟎
            𝟎 𝟎 𝟎 𝟎 c₄₄ 𝟎
            𝟎 𝟎 𝟎 𝟎 𝟎 c₄₄
        ],
    )
end
function (::ElasticConstantFitter{Cubic})(
    stresses::AbstractVector{<:EngineeringStress},
    strains::AbstractVector{<:EngineeringStrain},
)
    σ, ϵ = first(stresses), first(strains)
    σ₁, σ₂, σ₃ = σ[1:3]
    Aᵀ = [
        σ₁ σ₂ σ₃
        σ₂+σ₃ σ₁+σ₃ σ₂+σ₁
    ]
    s₁₁, s₁₂ = inv(Aᵀ * transpose(Aᵀ)) * Aᵀ * ϵ[1:3]  # If 𝐴 is well-conditioned, using the normal equations is around as accurate as other methods and is also the fastest. https://math.stackexchange.com/a/3252377/115512
    s₄₄ = dot(σ[4:6], ϵ[4:6]) / sum(abs2, σ[4:6])  # B = σ[4:6], s₄₄ = inv(Bᵀ * B) * Bᵀ * σ[4:6]
    𝟎 = zero(s₁₁)
    return ComplianceMatrix(
        [
            s₁₁ s₁₂ s₁₂ 𝟎 𝟎 𝟎
            s₁₂ s₁₁ s₁₂ 𝟎 𝟎 𝟎
            s₁₂ s₁₂ s₁₁ 𝟎 𝟎 𝟎
            𝟎 𝟎 𝟎 s₄₄ 𝟎 𝟎
            𝟎 𝟎 𝟎 𝟎 s₄₄ 𝟎
            𝟎 𝟎 𝟎 𝟎 𝟎 s₄₄
        ],
    )
end
function (x::ElasticConstantFitter)(
    strains::AbstractVector{<:TensorStrain},
    stresses::AbstractVector{<:TensorStress},
)
    c = x(EngineeringStrain.(strains), EngineeringStress.(stresses))
    return StiffnessTensor(c)
end
function (x::ElasticConstantFitter)(
    stresses::AbstractVector{<:TensorStress},
    strains::AbstractVector{<:TensorStrain},
)
    s = x(EngineeringStress.(stresses), EngineeringStrain.(strains))
    return ComplianceTensor(s)
end
for (X, Y) in ((:EngineeringStrain, :EngineeringStress), (:TensorStrain, :TensorStress))
    @eval begin
        (x::ElasticConstantFitter)(
            strains::AbstractVector{<:$X},
            stresses::AbstractVector{<:$Y},
            σ₀::$Y,
        ) = x(strains, map(Base.Fix2(-, σ₀), stresses))  # Subtract a common initial value σ₀
        (x::ElasticConstantFitter)(
            stresses::AbstractVector{<:$Y},
            strains::AbstractVector{<:$X},
            ϵ₀::$X,
        ) = x(stresses, map(Base.Fix2(-, ϵ₀), strains))
    end
end
