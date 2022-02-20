using Compat: only
using Crystallography: Triclinic, Trigonal

struct ElasticConstantSolver{T<:CrystalSystem}
    system::T
end

function (::ElasticConstantSolver{Triclinic})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    cᵢⱼ = [_cᵢⱼ(strains, stresses, i, j) for i in 1:6 for j in i:6]
    return StiffnessMatrix(cᵢⱼ...)
end
function (::ElasticConstantSolver{Orthorhombic})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    c₁₁ = _cᵢⱼ(strains, stresses, 1, 1)
    c₁₂ = _cᵢⱼ(strains, stresses, 1, 2)
    c₁₃ = _cᵢⱼ(strains, stresses, 1, 3)
    c₂₂ = _cᵢⱼ(strains, stresses, 2, 2)
    c₂₃ = _cᵢⱼ(strains, stresses, 2, 3)
    c₃₃ = _cᵢⱼ(strains, stresses, 3, 3)
    c₄₄ = _cᵢⱼ(strains, stresses, 4, 4)
    c₅₅ = _cᵢⱼ(strains, stresses, 5, 5)
    c₆₆ = _cᵢⱼ(strains, stresses, 6, 6)
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
function (::ElasticConstantSolver{Tetragonal})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    c₁₁ = _cᵢⱼ(strains, stresses, 1, 1)
    c₁₂ = _cᵢⱼ(strains, stresses, 1, 2)
    c₁₃ = _cᵢⱼ(strains, stresses, 1, 3)
    c₁₆ = _cᵢⱼ(strains, stresses, 1, 6)
    c₂₂ = _cᵢⱼ(strains, stresses, 2, 2)
    c₃₃ = _cᵢⱼ(strains, stresses, 3, 3)
    c₄₄ = _cᵢⱼ(strains, stresses, 4, 4)
    c₆₆ = _cᵢⱼ(strains, stresses, 6, 6)
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        [
            c₁₁ c₁₂ c₁₃ 𝟎 𝟎 c₁₆
            c₁₂ c₂₂ c₁₃ 𝟎 𝟎 -c₁₆
            c₁₃ c₁₃ c₃₃ 𝟎 𝟎 𝟎
            𝟎 𝟎 𝟎 c₄₄ 𝟎 𝟎
            𝟎 𝟎 𝟎 𝟎 c₄₄ 𝟎
            c₁₆ -c₁₆ 𝟎 𝟎 𝟎 c₆₆
        ],
    )
end
function (::ElasticConstantSolver{Trigonal})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    c₁₁ = _cᵢⱼ(strains, stresses, 1, 1)
    c₁₂ = _cᵢⱼ(strains, stresses, 1, 2)
    c₁₃ = _cᵢⱼ(strains, stresses, 1, 3)
    c₁₄ = _cᵢⱼ(strains, stresses, 1, 4)
    c₁₅ = _cᵢⱼ(strains, stresses, 1, 5)
    c₂₂ = _cᵢⱼ(strains, stresses, 2, 2)
    c₃₃ = _cᵢⱼ(strains, stresses, 3, 3)
    c₄₄ = _cᵢⱼ(strains, stresses, 4, 4)
    𝟎 = zero(c₁₁)
    return StiffnessMatrix(
        [
            c₁₁ c₁₂ c₁₃ c₁₄ c₁₅ 𝟎
            c₁₂ c₂₂ c₁₃ -c₁₄ -c₁₅ 𝟎
            c₁₃ c₁₃ c₃₃ 𝟎 𝟎 𝟎
            𝟎 𝟎 𝟎 c₄₄ 𝟎 -c₁₅
            𝟎 𝟎 𝟎 𝟎 c₄₄ c₁₄
            𝟎 𝟎 𝟎 𝟎 𝟎 (c₁₁-c₁₂)/2
        ],
    )
end
function (::ElasticConstantSolver{Hexagonal})(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)
    c₁₁ = _cᵢⱼ(strains, stresses, 1, 1)
    c₁₂ = _cᵢⱼ(strains, stresses, 1, 2)
    c₁₃ = _cᵢⱼ(strains, stresses, 1, 3)
    c₃₃ = _cᵢⱼ(strains, stresses, 3, 3)
    c₄₄ = _cᵢⱼ(strains, stresses, 4, 4)
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
    c₁₁ = _cᵢⱼ(strains, stresses, 1, 1)
    c₁₂ = _cᵢⱼ(strains, stresses, 1, 2)
    c₄₄ = _cᵢⱼ(strains, stresses, 4, 4)
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

function _find_nonzero_element(strain::EngineeringStrain)
    indices = findall(!iszero, strain)
    return only(indices)
end

function _pick_from(strains::AbstractVector{<:EngineeringStrain})
    indices = map(_find_nonzero_element, strains)
    function _at_index(desired_index)
        positions = findall(==(desired_index), indices)  # No duplicated directions allowed
        return only(positions)
    end
end

function _cᵢⱼ(
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
    j,
    i,
)
    position = _pick_from(strains)(j)
    σᵢ, ϵⱼ = stresses[position][i], strains[position][j]
    return -σᵢ / ϵⱼ
end
