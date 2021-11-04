using LinearAlgebra: eigvals, eigvecs
using Tensorial: stress_invariants, deviatoric_stress_invariants, vol, dev

export principal_values,
    principal_axes, principal_invariants, main_invariants, hydrostatic, deviatoric

principal_values(x::Union{TensorStress,TensorStrain}) = eigvals(x)

principal_axes(x::Union{TensorStress,TensorStrain}) = eigvecs(x)

principal_invariants(x::Union{TensorStress,TensorStrain}) = stress_invariants(x.data)

main_invariants(x::Union{TensorStress,TensorStrain}) = deviatoric_stress_invariants(x.data)

hydrostatic(x::Union{TensorStress,TensorStrain}) = vol(x.data)

deviatoric(x::Union{TensorStress,TensorStrain}) = dev(x.data)
