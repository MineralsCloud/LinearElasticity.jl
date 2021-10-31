module LinearElasticity

using Tensorial: SymmetricSecondOrderTensor, SymmetricFourthOrderTensor, Vec

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
    data::SymmetricSecondOrderTensor{3,T,6}
end
struct TensorStrain{T} <: Strain{T,2}
    data::SymmetricSecondOrderTensor{3,T,6}
end
struct StiffnessTensor{T} <: Stiffness{T,4}
    data::SymmetricFourthOrderTensor{3,T}
end
struct ComplianceTensor{T} <: Compliance{T,4}
    data::SymmetricFourthOrderTensor{3,T}
end
struct EngineeringStress{T} <: Stress{T,1}
    data::Vec{6,T}
end
struct EngineeringStrain{T} <: Strain{T,1}
    data::Vec{6,T}
end
struct StiffnessMatrix{T} <: Stiffness{T,2}
    data::SymmetricSecondOrderTensor{6,T,21}
end
struct ComplianceMatrix{T} <: Compliance{T,2}
    data::SymmetricSecondOrderTensor{6,T,21}
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
    @eval begin
        $T(m::AbstractMatrix) = $T(SymmetricSecondOrderTensor{3}(m))
        $T(data...) = $T(SymmetricSecondOrderTensor{3}(data...))
    end
end
for T in (:EngineeringStress, :EngineeringStrain)
    @eval begin
        $T(v::AbstractVector) = $T(Vec{6}(v))
        $T(data...) = $T(Vec{6}(data...))
    end
end
for T in (:StiffnessMatrix, :ComplianceMatrix)
    @eval begin
        $T(m::AbstractMatrix) = $T(SymmetricSecondOrderTensor{6}(m))
        $T(data...) = $T(SymmetricSecondOrderTensor{6}(data...))
    end
end

include("conversion.jl")
include("symmetry_criteria.jl")
include("invariants.jl")
include("StabilityCriteria.jl")
include("Moduli.jl")

end
