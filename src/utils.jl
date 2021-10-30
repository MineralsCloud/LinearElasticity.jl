using LinearAlgebra: tr, det, eigvals, eigvecs

export principal_values, principal_axes, principal_invariants, main_invariants, issystem

issystem(C::CrystalSystem, x::Union{MatrixStiffness,MatrixCompliance}) =
    all(symmetry_criteria(C, x))

principal_values(x::Union{Stress,Strain}) = eigvals(x)

principal_axes(x::Union{Stress,Strain}) = eigvecs(x)

principal_invariants(x::Union{Stress,Strain}, n::Int) = principal_invariants(x, Val(n))
principal_invariants(x::Union{Stress,Strain}, ::Val{1}) = tr(x)
function principal_invariants(x::Union{Stress,Strain}, ::Val{2})
    λ1, λ2, λ3 = principal_values(x)
    return λ1 * λ2 + λ1 * λ3 + λ2 * λ3
end # function principal_invariants
principal_invariants(x::Union{Stress,Strain}, ::Val{3}) = det(x)

main_invariants(x::Union{Stress,Strain}, n::Int) = main_invariants(x, Val(n))
main_invariants(x::Union{Stress,Strain}, ::Val{1}) = principal_invariants(x, 1)
main_invariants(x::Union{Stress,Strain}, ::Val{2}) =
    principal_invariants(x, 1)^2 - 2 * principal_invariants(x, 2)
main_invariants(x::Union{Stress,Strain}, ::Val{3}) =
    principal_invariants(x, 1)^3 +
    3 *
    (principal_invariants(x, 3) - principal_invariants(x, 1) * principal_invariants(x, 2))
