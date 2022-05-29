module LinearElasticity

using ConstructionBase: constructorof
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

Base.parent(A::Union{Stress,Strain,Stiffness,Compliance}) = A.data

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

# See https://github.com/JuliaLang/julia/blob/cb9acf5/base/arraymath.jl#L19-L26
for op in (:*, :/)
    @eval Base.$op(A::Union{Stress,Strain,Stiffness,Compliance}, B::Number) =
        constructorof(typeof(A))(Base.broadcast_preserving_zero_d($op, A, B))
end
Base.:*(B::Number, A::Union{Stress,Strain,Stiffness,Compliance}) = A * B

Base.:-(A::Union{Stress,Strain,Stiffness,Compliance}) =
    constructorof(typeof(A))(Base.broadcast_preserving_zero_d(-, A))

for op in (:+, :-)
    for T in (
        :TensorStress,
        :TensorStrain,
        :EngineeringStress,
        :EngineeringStrain,
        :StiffnessMatrix,
        :ComplianceMatrix,
        :StiffnessTensor,
        :ComplianceTensor,
    )
        @eval Base.$op(A::$T, B::$T) = $T(Base.broadcast_preserving_zero_d($op, A, B))
    end
end

include("conversion.jl")
include("invariants.jl")
include("SymmetryCriteria.jl")
# include("StabilityCriteria.jl")
include("Isotropic.jl")
include("distort.jl")
include("fitting.jl")
include("calc.jl")

end
