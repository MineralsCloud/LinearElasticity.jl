using Compat: only
using Crystallography: Triclinic

struct ElasticConstantSolver{T<:CrystalSystem}
    system::T
end

function (::ElasticConstantSolver{Triclinic})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    cᵢⱼ = [_calculate_cij(strains, stresses, i, j) for i in 1:6 for j in i:6]
    return StiffnessMatrix(cᵢⱼ)
end
function (::ElasticConstantSolver{Orthorhombic})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    c₁₁ = _calculate_cij(strains, stresses, 1, 1)
    c₁₂ = _calculate_cij(strains, stresses, 1, 2)
    c₁₃ = _calculate_cij(strains, stresses, 1, 3)
    c₂₂ = _calculate_cij(strains, stresses, 2, 2)
    c₂₃ = _calculate_cij(strains, stresses, 2, 3)
    c₃₃ = _calculate_cij(strains, stresses, 3, 3)
    c₄₄ = _calculate_cij(strains, stresses, 4, 4)
    c₅₅ = _calculate_cij(strains, stresses, 5, 5)
    c₆₆ = _calculate_cij(strains, stresses, 6, 6)
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
    )
end
function (::ElasticConstantSolver{Hexagonal})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    c₁₁ = _calculate_cij(strains, stresses, 1, 1)
    c₁₂ = _calculate_cij(strains, stresses, 1, 2)
    c₁₃ = _calculate_cij(strains, stresses, 1, 3)
    c₃₃ = _calculate_cij(strains, stresses, 3, 3)
    c₄₄ = _calculate_cij(strains, stresses, 4, 4)
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
function (::ElasticConstantSolver{Cubic})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    c₁₁ = _calculate_cij(strains, stresses, 1, 1)
    c₁₂ = _calculate_cij(strains, stresses, 1, 2)
    c₄₄ = _calculate_cij(strains, stresses, 4, 4)
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

function _find_nonzero_element(strain_or_stress::Union{EngineeringStress,EngineeringStrain})
    indices = findall(!iszero, strain_or_stress)
    return only(indices)
end

function _pick_nonzero(strains_or_stresses::AbstractVector)
    indices = map(_find_nonzero_element, strains_or_stresses)
    function _at_index(desired_index)
        positions = findall(==(desired_index), indices)  # No duplicated directions allowed
        position = only(positions)
        return strains_or_stresses[position]
    end
end

_isnegative(number) = number < zero(number)

_cij(ϵᵢ₊, ϵᵢ₋, σⱼ₊, σⱼ₋) = (σⱼ₊ - σⱼ₋) / (ϵᵢ₊ - ϵᵢ₋)

function _calculate_cij(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
    i,
    j,
)
    ϵᵢ₊, ϵᵢ₋ = _pick_nonzero(strains)(i)
    σⱼ₊, σⱼ₋ = _pick_nonzero(stresses)(j)
    return _cij(ϵᵢ₊, ϵᵢ₋, σⱼ₊, σⱼ₋)
end
