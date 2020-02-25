module Elasticity

using Tensors: SymmetricTensor

export TensorStress,
    TensorStrain,
    TensorStiffness,
    TensorCompliance,
    EngineeringStress,
    EngineeringStrain,
    EngineeringCompliance,
    EngineeringStiffness

struct TensorStress{T} <: AbstractMatrix{T}
    data::SymmetricTensor{2,3,T}
end

struct TensorStrain{T} <: AbstractMatrix{T}
    data::SymmetricTensor{2,3,T}
end

struct TensorStiffness{T} <: AbstractArray{T,4}
    data::SymmetricTensor{4,3,T}
end

struct TensorCompliance{T} <: AbstractArray{T,4}
    data::SymmetricTensor{4,3,T}
end

struct EngineeringStress{T} <: AbstractVector{T}
    data::NTuple{6,T}
end

struct EngineeringStrain{T} <: AbstractVector{T}
    data::NTuple{6,T}
end

struct EngineeringStiffness{T} <: AbstractMatrix{T}
    data::SymmetricTensor{2,6,T}
end

struct EngineeringCompliance{T} <: AbstractMatrix{T}
    data::SymmetricTensor{2,6,T}
end

const Stress = Union{TensorStress,EngineeringStress}
const Strain = Union{TensorStrain,EngineeringStrain}
const Stiffness = Union{TensorStiffness,EngineeringStiffness}
const Compliance = Union{TensorCompliance,EngineeringCompliance}

Base.size(::Union{TensorStress,TensorStrain}) = (3, 3)
Base.size(::Union{TensorStiffness,TensorCompliance}) = (3, 3, 3, 3)
Base.size(::Union{EngineeringStress,EngineeringStrain}) = (6,)
Base.size(::Union{EngineeringStiffness,EngineeringCompliance}) = (6, 6)

Base.getindex(A::Union{Stress,Strain,Stiffness,Compliance}, I::Vararg{Int}) =
    getindex(A.data, I...)

Base.inv(c::TensorStiffness) = TensorCompliance(inv(c))
Base.inv(s::TensorCompliance) = TensorStiffness(inv(s))
Base.inv(c::EngineeringStiffness) = EngineeringCompliance(inv(c))
Base.inv(s::EngineeringCompliance) = EngineeringStiffness(inv(s))

for T in (:TensorStress, :TensorStrain)
    eval(quote
        $T(m::AbstractMatrix) = $T(SymmetricTensor{2,3}(m))
        $T(v::Union{AbstractVector,NTuple{6}}) = $T(SymmetricTensor{2,3}(v))
        $T(xx, xy, xz, yy, yz, zz) = $T((xx, xy, xz, yy, yz, zz))
    end)
end

end
