module LinearElasticity

using LinearAlgebra: tr, det, eigvals, eigvecs, issymmetric

using Crystallography:
    CrystalSystem, Cubic, Hexagonal, Tetragonal, Trigonal, Orthorhombic, Monoclinic
using StaticArrays: SHermitianCompact, SArray, SMatrix, SVector

export TensorStress,
    TensorStrain,
    TensorStiffness,
    TensorCompliance,
    EngineeringStress,
    EngineeringStrain,
    EngineeringCompliance,
    EngineeringStiffness
export principal_values, principal_axes, principal_invariants, main_invariants, issystem

abstract type Stress{T,N} <: AbstractArray{T,N} end
abstract type Strain{T,N} <: AbstractArray{T,N} end
abstract type Stiffness{T,N} <: AbstractArray{T,N} end
abstract type Compliance{T,N} <: AbstractArray{T,N} end
struct TensorStress{T} <: Stress{T,2}
    data::SHermitianCompact{3,T}
end
struct TensorStrain{T} <: Strain{T,2}
    data::SHermitianCompact{3,T}
end
struct TensorStiffness{T} <: Stiffness{T,4}
    data::SArray{Tuple{3,3,3,3},T}
end
struct TensorCompliance{T} <: Compliance{T,4}
    data::SArray{Tuple{3,3,3,3},T}
end
struct EngineeringStress{T} <: Stress{T,1}
    data::SVector{6,T}
end
struct EngineeringStrain{T} <: Strain{T,1}
    data::SVector{6,T}
end
struct EngineeringStiffness{T} <: Stiffness{T,2}
    data::SHermitianCompact{6,T}
end
struct EngineeringCompliance{T} <: Compliance{T,2}
    data::SHermitianCompact{6,T}
end

issystem(C::CrystalSystem, x::Union{EngineeringStiffness,EngineeringCompliance}) =
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

Base.size(::Union{TensorStress,TensorStrain}) = (3, 3)
Base.size(::Union{TensorStiffness,TensorCompliance}) = (3, 3, 3, 3)
Base.size(::Union{EngineeringStress,EngineeringStrain}) = (6,)
Base.size(::Union{EngineeringStiffness,EngineeringCompliance}) = (6, 6)

Base.getindex(A::Union{EngineeringStress,EngineeringStrain}, i::Int) = getindex(A.data, i)
Base.getindex(
    A::Union{
        TensorStress,
        TensorStrain,
        TensorStiffness,
        TensorCompliance,
        EngineeringStiffness,
        EngineeringCompliance,
    },
    I::Vararg{Int},
) = getindex(A.data, I...)

for T in (:TensorStress, :TensorStrain)
    # See https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#StaticArrays.SHermitianCompact
    @eval begin
        $T(m::AbstractMatrix) = $T(SHermitianCompact{3}(m))
        $T(v::AbstractVector) = $T(SHermitianCompact(SVector{6}(v)))
        $T(t::NTuple{9}) = $T(SHermitianCompact{3}(t))
    end
end
for T in (:EngineeringStress, :EngineeringStrain)
    # See https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#StaticArrays.SHermitianCompact
    @eval begin
        $T(v::AbstractVector) = $T(SVector{6}(v))
    end
end
for T in (:EngineeringStiffness, :EngineeringCompliance)
    @eval begin
        $T(m::AbstractMatrix) = $T(SHermitianCompact{6}(m))
        $T(v::AbstractVector) = $T(SHermitianCompact(SVector{21}(v)))
        $T(t::NTuple{36}) = $T(SHermitianCompact{6}(t))
    end
end
for T in (:TensorStiffness, :TensorCompliance)
    @eval begin
        $T(a::AbstractArray) = $T(SArray{Tuple{3,3,3,3}}(a))
    end
end

const VOIGT_INDICES =
    ((1, 1), (2, 2), (3, 3), (3, 2), (3, 1), (2, 1), (2, 3), (1, 3), (1, 2))

include("conversion.jl")
include("symmetry_criteria.jl")
include("StabilityCriteria.jl")
include("Moduli.jl")

end
