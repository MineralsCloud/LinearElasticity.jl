struct ElasticConstantSolver{T<:CrystalSystem}
    system::T
end

function (::ElasticConstantSolver{Cubic})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    indices = map(_whichindex, strains)
    strain1 = _select(strains, indices, 1)
    stress1 = _select(stresses, indices, 1)
    c₁₁ = _cij(strain1[1], strain1[2], stress1[1], stress1[2])
    stress2 = _select(stresses, indices, 2)
    c₁₂ = _cij(strain1[1], strain1[2], stress2[1], stress2[2])
    strain4 = _select(strains, indices, 4)
    stress4 = _select(stresses, indices, 4)
    c₄₄ = _cij(strain4[1], strain4[2], stress4[1], stress4[2])
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
