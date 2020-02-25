module Elasticity

using Tensors: SymmetricTensor

export TensorStress, TensorStrain

struct TensorStress{T}
    data::SymmetricTensor{2,3,T}
end

struct TensorStrain{T}
    data::SymmetricTensor{2,3,T}
end

for T in (:TensorStress, :TensorStrain)
    eval(quote
        $T(v::Union{AbstractVector,NTuple{6}}) = $T(SymmetricTensor{2,3}(v))
        $T(xx, xy, xz, yy, yz, zz) = $T((xx, xy, xz, yy, yz, zz))
    end)
end

end
