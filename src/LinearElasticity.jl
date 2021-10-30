module LinearElasticity

using StaticArrays: SHermitianCompact, SMatrix, SVector
using Tensorial: SymmetricFourthOrderTensor

export TensorStress,
    TensorStrain,
    StiffnessTensor,
    ComplianceTensor,
    EngineeringStress,
    EngineeringStrain,
    ComplianceMatrix,
    StiffnessMatrix

abstract type Stress{T,N} <: AbstractArray{T,N} end
abstract type Strain{T,N} <: AbstractArray{T,N} end
abstract type Stiffness{T,N} <: AbstractArray{T,N} end
abstract type Compliance{T,N} <: AbstractArray{T,N} end
struct TensorStress{T} <: Stress{T,2}
    data::SHermitianCompact{3,T}
end
struct TensorStrain{T} <: Strain{T,2}
    data::SHermitianCompact{3,T}
end
struct StiffnessTensor{T} <: Stiffness{T,4}
    data::SymmetricFourthOrderTensor{3,T}
end
struct ComplianceTensor{T} <: Compliance{T,4}
    data::SymmetricFourthOrderTensor{3,T}
end
struct EngineeringStress{T} <: Stress{T,1}
    data::SVector{6,T}
end
struct EngineeringStrain{T} <: Strain{T,1}
    data::SVector{6,T}
end
struct StiffnessMatrix{T} <: Stiffness{T,2}
    data::SHermitianCompact{6,T}
end
struct ComplianceMatrix{T} <: Compliance{T,2}
    data::SHermitianCompact{6,T}
end

TensorStress(x::EngineeringStress) = convert(TensorStress{eltype(x)}, x)

TensorStrain(x::EngineeringStrain) = convert(TensorStrain{eltype(x)}, x)

StiffnessTensor(x::StiffnessMatrix) = convert(StiffnessTensor{eltype(x)}, x)

ComplianceTensor(x::ComplianceMatrix) = convert(ComplianceTensor{eltype(x)}, x)

EngineeringStress(x::TensorStress) = convert(EngineeringStress{eltype(x)}, x)

EngineeringStrain(x::TensorStrain) = convert(EngineeringStrain{eltype(x)}, x)

StiffnessMatrix(x::StiffnessTensor) = convert(StiffnessMatrix{eltype(x)}, x)

ComplianceMatrix(x::ComplianceTensor) = convert(ComplianceMatrix{eltype(x)}, x)

Base.size(::Union{TensorStress,TensorStrain}) = (3, 3)
Base.size(::Union{StiffnessTensor,ComplianceTensor}) = (3, 3, 3, 3)
Base.size(::Union{EngineeringStress,EngineeringStrain}) = (6,)
Base.size(::Union{StiffnessMatrix,ComplianceMatrix}) = (6, 6)

Base.getindex(A::Union{Stress,Strain,Stiffness,Compliance}, i...) = getindex(A.data, i...)

Base.IndexStyle(::Type{<:Union{Stress,Strain,Stiffness,Compliance}}) = IndexLinear()

for T in (:TensorStress, :TensorStrain)
    # See https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#StaticArrays.SHermitianCompact
    @eval begin
        $T(m::AbstractMatrix) = $T(SHermitianCompact{3}(m))
        $T(v::AbstractVector) = $T(SHermitianCompact(SVector{6}(v)))
        $T(t::NTuple{9}) = $T(SHermitianCompact{3}(t))
    end
end
for T in (:EngineeringStress, :EngineeringStrain)
    # See https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#StaticArrays.SHermitianCompact
    @eval begin
        $T(v::AbstractVector) = $T(SVector{6}(v))
    end
end
for T in (:StiffnessMatrix, :ComplianceMatrix)
    @eval begin
        $T(m::AbstractMatrix) = $T(SHermitianCompact{6}(m))
        $T(v::AbstractVector) = $T(SHermitianCompact(SVector{21}(v)))
        $T(t::NTuple{36}) = $T(SHermitianCompact{6}(t))
    end
end

include("conversion.jl")
include("symmetry_criteria.jl")
include("utils.jl")
include("StabilityCriteria.jl")
include("Moduli.jl")

end
