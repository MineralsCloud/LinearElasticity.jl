module Elasticity

using Tensors: SymmetricTensor

export TensorStress, TensorStrain, TensorStiffness, TensorCompliance

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

Base.size(::Union{TensorStress,TensorStrain}) = (3, 3)
Base.size(::Union{TensorStiffness,TensorCompliance}) = (3, 3, 3, 3)

Base.getindex(
    A::Union{TensorStress,TensorStrain,TensorStiffness,TensorCompliance},
    I::Vararg{Int},
) = getindex(A.data, I...)

for T in (:TensorStress, :TensorStrain)
    eval(quote
        $T(m::AbstractMatrix) = $T(SymmetricTensor{2,3}(m))
        $T(v::Union{AbstractVector,NTuple{6}}) = $T(SymmetricTensor{2,3}(v))
        $T(xx, xy, xz, yy, yz, zz) = $T((xx, xy, xz, yy, yz, zz))
    end)
end

end
