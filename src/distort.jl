using CrystallographyBase: Lattice, CrystalSystem, Cubic, Hexagonal
using LinearAlgebra: I, norm, dot

export ElasticConstantFitter, distort

# See https://link.springer.com/content/pdf/10.1007%2F978-3-7091-0382-1_7.pdf and https://doi.org/10.2138/am-1997-1-207
distort(lattice::Lattice, strain::TensorStrain) = Lattice((I + strain.data) * lattice.data)
distort(lattice::Lattice, strain::EngineeringStrain) =
    distort(lattice, TensorStrain(strain))

struct ElasticConstantFitter{T<:CrystalSystem}
    system::T
end

function (::ElasticConstantFitter{Cubic})(ϵ::EngineeringStrain, σ::EngineeringStress)
    σ = map(Base.Fix1(oftype, σ[1]) ∘ float, σ)
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
function (x::ElasticConstantFitter)(ϵ::TensorStrain, σ::TensorStress)
    c = x(EngineeringStrain(ϵ), EngineeringStress(σ))
    return StiffnessTensor(c)
end
