using CrystallographyBase: Lattice, CrystalSystem, Cubic
using LinearAlgebra: I, svd, diagm, qr

export distort, fit

# See https://link.springer.com/content/pdf/10.1007%2F978-3-7091-0382-1_7.pdf and https://doi.org/10.2138/am-1997-1-207
distort(lattice::Lattice, strain::TensorStrain) = Lattice((I + strain.data) * lattice.data)
distort(lattice::Lattice, strain::EngineeringStrain) =
    distort(lattice, TensorStrain(strain))

function fit(ϵ::EngineeringStrain, σ::EngineeringStress, ::Cubic)
    σ = map(Base.Fix1(oftype, σ[1]) ∘ float, σ)
    ϵ₁, ϵ₂, ϵ₃ = ϵ[1:3]
    A = [
        ϵ₁ ϵ₂+ϵ₃
        ϵ₂ ϵ₁+ϵ₃
        ϵ₃ ϵ₂+ϵ₁
    ]
    # https://discourse.julialang.org/t/why-does-julia-systematically-fails-when-doing-operation/67242/9
    Aᵀ = transpose(A)
    c11, c12 = inv(Aᵀ * A) * Aᵀ * σ[1:3]  # If 𝐴 is well-conditioned, using the normal equations is around as accurate as other methods and is also the fastest. https://math.stackexchange.com/a/3252377/115512
    c44 = first(curve_fit(_f, ϵ[4:6], σ[4:6], [σ[4] / ϵ[4]]).param)
    𝟘 = zero(c11)
    data =
        [
            c11 c12 c12 𝟘 𝟘 𝟘
            c12 c11 c12 𝟘 𝟘 𝟘
            c12 c12 c11 𝟘 𝟘 𝟘
            𝟘 𝟘 𝟘 c44 𝟘 𝟘
            𝟘 𝟘 𝟘 𝟘 c44 𝟘
            𝟘 𝟘 𝟘 𝟘 𝟘 c44
        ] * oneunit(σ[1])
    return StiffnessMatrix(data)
end
function fit(ϵ::TensorStrain, σ::TensorStress, x::CrystalSystem)
    c = fit(EngineeringStrain(ϵ), EngineeringStress(σ), x)
    return StiffnessTensor(c)
end

_f(x, p) = p[1] .* x
