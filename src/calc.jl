struct ElasticConstantSolver{T<:CrystalSystem}
    system::T
end

function (::ElasticConstantSolver{Cubic})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    indices = map(_whichindex, strains)
    ϵ₁₊, ϵ₁₋ = _select(strains, indices, 1)
    σ₁₊, σ₁₋ = _select(stresses, indices, 1)
    c₁₁ = _cij(ϵ₁₊, ϵ₁₋, σ₁₊, σ₁₋)
    σ₂₊, σ₂₋ = _select(stresses, indices, 2)
    c₁₂ = _cij(ϵ₁₊, ϵ₁₋, σ₂₊, σ₂₋)
    ϵ₄₊, ϵ₄₋ = _select(strains, indices, 4)
    σ₄₊, σ₄₋ = _select(stresses, indices, 4)
    c₄₄ = _cij(ϵ₄₊, ϵ₄₋, σ₄₊, σ₄₋)
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

_select(v, indices, index) = v[filter(==(index), indices)]

_cij(ϵᵢ₊, ϵᵢ₋, σⱼ₊, σⱼ₋) = (σⱼ₊ - σⱼ₋) / (ϵᵢ₊ - ϵᵢ₋)

function _whichindex(x)
    indices = findall(!iszero, x)
    @assert length(indices) == 1
    return first(indices)
end
