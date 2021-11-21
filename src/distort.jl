using CrystallographyBase: Lattice, CrystalSystem, Cubic, Hexagonal
using LinearAlgebra: I, norm, dot

export ElasticConstantFitter, distortby, distort, strainstate

# See https://link.springer.com/content/pdf/10.1007%2F978-3-7091-0382-1_7.pdf and https://doi.org/10.2138/am-1997-1-207
distortby(lattice::Lattice, strain::TensorStrain) =
    Lattice((I + strain.data) * lattice.data)
distortby(lattice::Lattice, strain::EngineeringStrain) =
    distortby(lattice, TensorStrain(strain))
const distort = distortby  # For the sake of compatibility

strainstate(old::Lattice, new::Lattice) = TensorStrain(new.data / old.data - I)

struct ElasticConstantFitter{T<:CrystalSystem}
    system::T
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
    𝟘 = zero(c₁₁)
    return StiffnessMatrix(
        [
            c₁₁ c₁₂ c₁₃ 𝟘 𝟘 𝟘
            c₁₂ c₁₁ c₁₃ 𝟘 𝟘 𝟘
            c₁₃ c₁₃ c₃₃ 𝟘 𝟘 𝟘
            𝟘 𝟘 𝟘 c₄₄ 𝟘 𝟘
            𝟘 𝟘 𝟘 𝟘 c₄₄ 𝟘
            𝟘 𝟘 𝟘 𝟘 𝟘 c₆₆
        ],
    )
end
function (::ElasticConstantFitter{Cubic})(ϵ::EngineeringStrain, σ::EngineeringStress)
    ϵ₁, ϵ₂, ϵ₃ = ϵ[1:3]
    Aᵀ = [
        ϵ₁ ϵ₂ ϵ₃
        ϵ₂+ϵ₃ ϵ₁+ϵ₃ ϵ₂+ϵ₁
    ]
    c₁₁, c₁₂ = inv(Aᵀ * transpose(Aᵀ)) * Aᵀ * σ[1:3]  # If 𝐴 is well-conditioned, using the normal equations is around as accurate as other methods and is also the fastest. https://math.stackexchange.com/a/3252377/115512
    c₄₄ = dot(ϵ[4:6], σ[4:6]) / sum(abs2, ϵ[4:6])  # B = ϵ[4:6], c₄₄ = inv(Bᵀ * B) * Bᵀ * σ[4:6]
    𝟘 = zero(c₁₁)
    return StiffnessMatrix(
        [
            c₁₁ c₁₂ c₁₂ 𝟘 𝟘 𝟘
            c₁₂ c₁₁ c₁₂ 𝟘 𝟘 𝟘
            c₁₂ c₁₂ c₁₁ 𝟘 𝟘 𝟘
            𝟘 𝟘 𝟘 c₄₄ 𝟘 𝟘
            𝟘 𝟘 𝟘 𝟘 c₄₄ 𝟘
            𝟘 𝟘 𝟘 𝟘 𝟘 c₄₄
        ],
    )
end
function (::ElasticConstantFitter{Cubic})(σ::EngineeringStress, ϵ::EngineeringStrain)
    σ₁, σ₂, σ₃ = σ[1:3]
    Aᵀ = [
        σ₁ σ₂ σ₃
        σ₂+σ₃ σ₁+σ₃ σ₂+σ₁
    ]
    s₁₁, s₁₂ = inv(Aᵀ * transpose(Aᵀ)) * Aᵀ * ϵ[1:3]  # If 𝐴 is well-conditioned, using the normal equations is around as accurate as other methods and is also the fastest. https://math.stackexchange.com/a/3252377/115512
    s₄₄ = dot(σ[4:6], ϵ[4:6]) / sum(abs2, σ[4:6])  # B = σ[4:6], s₄₄ = inv(Bᵀ * B) * Bᵀ * σ[4:6]
    𝟘 = zero(s₁₁)
    return ComplianceMatrix(
        [
            s₁₁ s₁₂ s₁₂ 𝟘 𝟘 𝟘
            s₁₂ s₁₁ s₁₂ 𝟘 𝟘 𝟘
            s₁₂ s₁₂ s₁₁ 𝟘 𝟘 𝟘
            𝟘 𝟘 𝟘 s₄₄ 𝟘 𝟘
            𝟘 𝟘 𝟘 𝟘 s₄₄ 𝟘
            𝟘 𝟘 𝟘 𝟘 𝟘 s₄₄
        ],
    )
end
function (x::ElasticConstantFitter)(ϵ::TensorStrain, σ::TensorStress)
    c = x(EngineeringStrain(ϵ), EngineeringStress(σ))
    return StiffnessTensor(c)
end
function (x::ElasticConstantFitter)(σ::TensorStress, ϵ::TensorStrain)
    s = x(EngineeringStress(σ), EngineeringStrain(ϵ))
    return ComplianceTensor(s)
end
for (X, Y) in ((:EngineeringStrain, :EngineeringStress), (:TensorStrain, :TensorStress))
    @eval begin
        (x::ElasticConstantFitter)(ϵ::$X, σ::$Y, σ₀::$Y) = x(ϵ, σ - σ₀)
        (x::ElasticConstantFitter)(σ::$Y, ϵ::$X, ϵ₀::$X) = x(σ, ϵ - ϵ₀)
    end
end
