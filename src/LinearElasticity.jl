module LinearElasticity

using LinearAlgebra: tr, det, eigvals, eigvecs

using StaticArrays: SHermitianCompact, SArray, SMatrix, SVector

export TensorStress,
    TensorStrain,
    TensorStiffness,
    TensorCompliance,
    EngineeringStress,
    EngineeringStrain,
    EngineeringCompliance,
    EngineeringStiffness
export principal_values, principal_axes, principal_invariants, main_invariants

struct TensorStress{T} <: AbstractMatrix{T}
    data::SHermitianCompact{3,T}
end

struct TensorStrain{T} <: AbstractMatrix{T}
    data::SHermitianCompact{3,T}
end

struct TensorStiffness{T} <: AbstractArray{T,4}
    data::SArray{Tuple{3,3,3,3},T}
end

struct TensorCompliance{T} <: AbstractArray{T,4}
    data::SArray{Tuple{3,3,3,3},T}
end

struct EngineeringStress{T} <: AbstractVector{T}
    data::SVector{6,T}
end

struct EngineeringStrain{T} <: AbstractVector{T}
    data::SVector{6,T}
end

struct EngineeringStiffness{T} <: AbstractMatrix{T}
    data::SHermitianCompact{6,T}
end

struct EngineeringCompliance{T} <: AbstractMatrix{T}
    data::SHermitianCompact{6,T}
end

const Stress = Union{TensorStress,EngineeringStress}
const Strain = Union{TensorStrain,EngineeringStrain}
const Stiffness = Union{TensorStiffness,EngineeringStiffness}
const Compliance = Union{TensorCompliance,EngineeringCompliance}

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

Base.inv(c::TensorStiffness) = TensorCompliance(inv(c))
Base.inv(s::TensorCompliance) = TensorStiffness(inv(s))
Base.inv(c::EngineeringStiffness) = EngineeringCompliance(inv(c))
Base.inv(s::EngineeringCompliance) = EngineeringStiffness(inv(s))

function Base.convert(::Type{TensorStress{T}}, s::EngineeringStress{T}) where {T}
    return TensorStress([s[1], s[6], s[5], s[2], s[4], s[3]])
end # function Base.convert
function Base.convert(::Type{EngineeringStress{T}}, s::TensorStress{T}) where {T}
    return EngineeringStress([s[1, 1], s[2, 2], s[3, 3], s[2, 3], s[1, 3], s[1, 2]])
end # function Base.convert
function Base.convert(::Type{TensorStrain{T}}, e::EngineeringStrain{T}) where {T}
    return TensorStrain([e[1], e[6] / 2, e[5] / 2, e[2], e[4] / 2, e[3]])
end # function Base.convert
function Base.convert(::Type{EngineeringStrain{T}}, e::TensorStrain{T}) where {T}
    return EngineeringStrain([
        e[1, 1],
        e[2, 2],
        e[3, 3],
        2 * e[2, 3],
        2 * e[1, 3],
        2 * e[1, 2],
    ])
end # function Base.convert
function Base.convert(::Type{EngineeringStiffness{T}}, c::TensorStiffness{T}) where {T}
    p = pairs(VOIGT_INDICES)
    return EngineeringStiffness(SymmetricTensor{2,6}((i, j) -> (println([p[i]..., p[j]...]); c[p[i]..., p[j]...])))
end # function Base.convert
function Base.convert(::Type{TensorStiffness{T}}, c::EngineeringStiffness{T}) where {T}
    d = Dict(zip(VOIGT_INDICES, 1:6))
    return TensorStiffness(SymmetricTensor{4,3}((i, j, k, l) -> c[d[(i, j)], d[(k, l)]]))
end # function Base.convert

for T in (:TensorStress, :TensorStrain)
    # See https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#StaticArrays.SHermitianCompact
    eval(quote
        $T(m::AbstractMatrix) = $T(SHermitianCompact{3}(m))
        $T(v::AbstractVector) = $T(SHermitianCompact(SVector{6}(v)))
        $T(t::NTuple{9}) = $T(SHermitianCompact{3}(t))
    end)
end
for T in (:EngineeringStress, :EngineeringStrain)
    # See https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#StaticArrays.SHermitianCompact
    eval(quote
        $T(v::AbstractVector) = $T(SVector{6}(v))
    end)
end
for T in (:EngineeringStiffness, :EngineeringCompliance)
    eval(quote
        $T(m::AbstractMatrix) = $T(SHermitianCompact{6}(m))
        $T(v::AbstractVector) = $T(SHermitianCompact(SVector{21}(v)))
        $T(t::NTuple{36}) = $T(SHermitianCompact{6}(t))
    end)
end

const VOIGT_INDICES = ((1, 1), (2, 2), (3, 3), (3, 2), (3, 1), (2, 1))

include("StabilityConditions.jl")

end
