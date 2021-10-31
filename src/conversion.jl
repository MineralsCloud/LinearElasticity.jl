using Tensorial: fromvoigt, ⊡

import Tensorial: contraction, double_contraction

export contraction, double_contraction

Base.inv(c::StiffnessTensor) = ComplianceTensor(inv(c.data))
Base.inv(s::ComplianceTensor) = StiffnessTensor(inv(s.data))
Base.inv(c::StiffnessMatrix) = ComplianceMatrix(inv(c.data))
Base.inv(s::ComplianceMatrix) = StiffnessMatrix(inv(s.data))

Base.convert(::Type{TensorStress{T}}, s::EngineeringStress{T}) where {T} =
    TensorStress(s[1], s[6], s[5], s[2], s[4], s[3])
Base.convert(::Type{EngineeringStress{T}}, σ::TensorStress{T}) where {T} =
    EngineeringStress(σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2])
Base.convert(::Type{TensorStrain{T}}, ϵ::EngineeringStrain{T}) where {T} =
    TensorStrain(ϵ[1], ϵ[6] / 2, ϵ[5] / 2, ϵ[2], ϵ[4] / 2, ϵ[3])
Base.convert(::Type{EngineeringStrain{T}}, ε::TensorStrain{T}) where {T} =
    EngineeringStrain(ε[1, 1], ε[2, 2], ε[3, 3], 2ε[2, 3], 2ε[1, 3], 2ε[1, 2])
Base.convert(::Type{StiffnessMatrix{T}}, c::StiffnessTensor{T}) where {T} = StiffnessMatrix(
    c[1, 1, 1, 1],
    c[1, 1, 2, 2],
    c[1, 1, 3, 3],
    c[1, 1, 2, 3],
    c[1, 1, 1, 3],
    c[1, 1, 1, 2],
    c[2, 2, 2, 2],
    c[2, 2, 3, 3],
    c[2, 2, 2, 3],
    c[2, 2, 1, 3],
    c[2, 2, 1, 2],
    c[3, 3, 3, 3],
    c[3, 3, 2, 3],
    c[3, 3, 1, 3],
    c[3, 3, 1, 2],
    c[2, 3, 2, 3],
    c[2, 3, 1, 3],
    c[2, 3, 1, 2],
    c[1, 3, 1, 3],
    c[1, 3, 1, 2],
    c[1, 2, 1, 2],
)
Base.convert(::Type{ComplianceMatrix{T}}, s::ComplianceTensor{T}) where {T} =
    ComplianceMatrix(
        s[1, 1, 1, 1],
        s[1, 1, 2, 2],
        s[1, 1, 3, 3],
        2s[1, 1, 2, 3],
        2s[1, 1, 1, 3],
        2s[1, 1, 1, 2],
        s[2, 2, 2, 2],
        s[2, 2, 3, 3],
        2s[2, 2, 2, 3],
        2s[2, 2, 1, 3],
        2s[2, 2, 1, 2],
        s[3, 3, 3, 3],
        2s[3, 3, 2, 3],
        2s[3, 3, 1, 3],
        2s[3, 3, 1, 2],
        4s[2, 3, 2, 3],
        4s[2, 3, 1, 3],
        4s[2, 3, 1, 2],
        4s[1, 3, 1, 3],
        4s[1, 3, 1, 2],
        4s[1, 2, 1, 2],
    )
Base.convert(::Type{StiffnessTensor{T}}, c::StiffnessMatrix{T}) where {T} =
    StiffnessTensor(fromvoigt(SymmetricFourthOrderTensor{3,T}, c.data))
function Base.convert(::Type{ComplianceTensor{T}}, s::ComplianceMatrix{T}) where {T}
    ComplianceTensor(
        SymmetricFourthOrderTensor{3,T}(function (i, j, k, l)
            if i == j && k == l
                return s[i, k]
            elseif i != j && k != l  # 4 = 9 - (2+3), 5 = 9 - (1+3), 6 = 9 - (1+2)
                return s[9-(i+j), 9-(k+l)] / 4
            elseif i == j && k != l
                return s[i, 9-(k+l)] / 2
            else  # i != j && k == l
                return s[9-(i+j), k] / 2
            end
        end),
    )
end

Base.:*(c::StiffnessMatrix, ϵ::EngineeringStrain) = EngineeringStress(c.data * ϵ.data)
Base.:*(s::ComplianceMatrix, σ::EngineeringStress) = EngineeringStrain(s.data * σ.data)

contraction(c::StiffnessTensor, ε::TensorStrain, ::Val{2}) = TensorStress(c.data ⊡ ε.data)
contraction(s::ComplianceTensor, σ::TensorStress, ::Val{2}) = TensorStrain(s.data ⊡ σ.data)

@inline double_contraction(c::StiffnessTensor, ε::TensorStrain) = contraction(c, ε, Val(2))
@inline double_contraction(s::ComplianceTensor, σ::TensorStress) = contraction(s, σ, Val(2))
