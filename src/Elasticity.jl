module Elasticity

using Tensors: SymmetricTensor, eigvals, eigvecs

export TensorStress,
    TensorStrain,
    TensorStiffness,
    TensorCompliance,
    EngineeringStress,
    EngineeringStrain,
    EngineeringCompliance,
    EngineeringStiffness
export principal_values, principal_axes

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

principal_values(x::Union{Stress,Strain}) = eigvals(x)
principal_axes(x::Union{Stress,Strain}) = eigvecs(x)

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

function Base.convert(::Type{<:TensorStress}, s::EngineeringStress)
    return TensorStress((s[1], s[6], s[5], s[2], s[4], s[3]))
end # function Base.convert
function Base.convert(::Type{<:EngineeringStress}, s::TensorStress)
    return EngineeringStress(s[1, 1], s[2, 2], s[3, 3], s[2, 3], s[1, 3], s[1, 2])
end # function Base.convert
function Base.convert(::Type{<:TensorStrain}, e::EngineeringStrain)
    return TensorStrain(e[1], e[6] / 2, e[5] / 2, e[2], e[4] / 2, e[3])
end # function Base.convert
function Base.convert(::Type{<:EngineeringStrain}, e::TensorStrain)
    return EngineeringStrain(
        e[1, 1],
        e[2, 2],
        e[3, 3],
        2 * e[2, 3],
        2 * e[1, 3],
        2 * e[1, 2],
    )
end # function Base.convert

for T in (:TensorStress, :TensorStrain)
    eval(quote
        $T(m::AbstractMatrix) = $T(SymmetricTensor{2,3}(m))
        $T(v::Union{AbstractVector,NTuple{6}}) = $T(SymmetricTensor{2,3}(v))
        $T(xx, xy, xz, yy, yz, zz) = $T((xx, xy, xz, yy, yz, zz))
    end)
end
for T in (:EngineeringStress, :EngineeringStrain)
    eval(quote
        $T(v::AbstractVector) = $T(Tuple(v))
        $T(xx, xy, xz, yy, yz, zz) = $T((xx, xy, xz, yy, yz, zz))
    end)
end
for T in (:EngineeringStiffness, :EngineeringCompliance)
    eval(quote
        $T(m::AbstractMatrix) = $T(SymmetricTensor{2,6}(m))
        $T(v::Union{AbstractVector,NTuple{21}}) = $T(SymmetricTensor{2,6}(v))
    end)
end

end
