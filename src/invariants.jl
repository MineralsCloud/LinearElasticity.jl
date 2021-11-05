using LinearAlgebra: Eigen, eigen
using Tensorial: stress_invariants, deviatoric_stress_invariants, vol, dev

export principal_values,
    principal_axes, principal_invariants, main_invariants, hydrostatic, deviatoric

principal_values(x::Union{TensorStress,TensorStrain}) = _eigen(x).values

principal_axes(x::Union{TensorStress,TensorStrain}) = _eigen(x).vectors

principal_invariants(x::Union{TensorStress,TensorStrain}) = stress_invariants(x.data)

main_invariants(x::Union{TensorStress,TensorStrain}) = deviatoric_stress_invariants(x.data)

for T in (:TensorStress, :TensorStrain)
    @eval begin
        hydrostatic(x::$T) = $T(vol(SymmetricSecondOrderTensor{3}(float(x))))
        deviatoric(x::$T) = $T(dev(SymmetricSecondOrderTensor{3}(float(x))))
    end
end

function _eigen(x)
    unitless = @. x / oneunit(x)
    eg = eigen(unitless)
    return Eigen(eg.values * oneunit(x[1]), eg.vectors)
end
