using LinearAlgebra: eigvals, eigvecs
using Tensorial: stress_invariants, deviatoric_stress_invariants

export principal_values, principal_axes, principal_invariants, main_invariants, issystem

issystem(C::CrystalSystem, x::Union{StiffnessMatrix,ComplianceMatrix}) =
    all(symmetry_criteria(C, x))

principal_values(x::Union{Stress,Strain}) = eigvals(x)

principal_axes(x::Union{Stress,Strain}) = eigvecs(x)

principal_invariants(x::Union{TensorStress,TensorStrain}) = stress_invariants(x.data)

main_invariants(x::Union{Stress,Strain}, n::Int) = main_invariants(x, Val(n))
main_invariants(x::Union{Stress,Strain}, ::Val{1}) = principal_invariants(x, 1)
main_invariants(x::Union{Stress,Strain}, ::Val{2}) =
    principal_invariants(x, 1)^2 - 2 * principal_invariants(x, 2)
main_invariants(x::Union{Stress,Strain}, ::Val{3}) =
    principal_invariants(x, 1)^3 +
    3 *
    (principal_invariants(x, 3) - principal_invariants(x, 1) * principal_invariants(x, 2))
