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

function symmetry_criteria(::Cubic, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        c[1, 1] == c[2, 2] == c[3, 3],
        c[4, 4] == c[5, 5] == c[6, 6],
        c[1, 2] == c[1, 3] == c[2, 3],
    ))
end
function symmetry_criteria(::Hexagonal, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        2c[6, 6] == c[1, 1] - c[1, 2],
    ))
end
function symmetry_criteria(::Tetragonal, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        if iszero(c[1, 6])  # Tetragonal (I) class, 4mm, -42m, 422, 4/mmm
            true
        else  # Tetragonal (II) class, 4, -4, 4/m
            c[1, 6] == -c[2, 6]
        end,
    ))
end
function symmetry_criteria(::Trigonal, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        2c[6, 6] == c[1, 1] - c[1, 2],
        c[1, 4] == -c[2, 4] == -c[5, 6],
        if iszero(c[1, 5])  # # Rhombohedral (I) class, 32, -3m, 3m
            true
        else  # Rhombohedral (II) class, 3, -3
            -c[1, 5] == c[2, 5] == c[4, 6]
        end,
    ))
end
function symmetry_criteria(::Orthorhombic, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        all(
            iszero,
            (
                c[1, 4],
                c[2, 4],
                c[3, 4],
                c[1, 5],
                c[2, 5],
                c[3, 5],
                c[4, 5],
                c[1, 6],
                c[2, 6],
                c[3, 6],
                c[4, 6],
                c[5, 6],
            ),
        ),
        all(
            !iszero,
            (
                c[1, 1],
                c[1, 2],
                c[1, 3],
                c[2, 2],
                c[2, 3],
                c[3, 3],
                c[4, 4],
                c[5, 5],
                c[6, 6],
            ),
        ),
    ))
end
function symmetry_criteria(::Monoclinic, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        all(
            iszero,
            (c[1, 4], c[2, 4], c[3, 4], c[4, 5], c[1, 6], c[2, 6], c[3, 6], c[5, 6]),
        ),
    ))
end
symmetry_criteria(C::CrystalSystem, s::EngineeringCompliance) = symmetry_criteria(C, inv(s))

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
    p, dim = pairs(VOIGT_INDICES), 6
    # From https://github.com/KristofferC/Tensors.jl/blob/bff451c/src/utilities.jl#L5-L14
    return EngineeringStiffness([c[p[i]..., p[j]...] for i = 1:dim, j = 1:dim])
end # function Base.convert
function Base.convert(::Type{TensorStiffness{T}}, c::EngineeringStiffness{T}) where {T}
    d, dim = Dict(zip(VOIGT_INDICES, [1, 2, 3, 4, 5, 6, 4, 5, 6])), 3
    return TensorStiffness([
        c[d[(i, j)], d[(k, l)]] for i = 1:dim, j = 1:dim, k = 1:dim, l = 1:dim
    ])
end # function Base.convert
function Base.convert(::Type{EngineeringCompliance{T}}, c::TensorCompliance{T}) where {T}  # FIXME Rules are wrong
    p, dim = pairs(VOIGT_INDICES), 6
    # From https://github.com/KristofferC/Tensors.jl/blob/bff451c/src/utilities.jl#L5-L14
    return EngineeringStiffness([c[p[i]..., p[j]...] for i = 1:dim, j = 1:dim])
end # function Base.convert
function Base.convert(::Type{TensorCompliance{T}}, c::EngineeringCompliance{T}) where {T}  # FIXME Rules are wrong
    d, dim = Dict(zip(VOIGT_INDICES, [1, 2, 3, 4, 5, 6, 4, 5, 6])), 3
    return TensorStiffness([
        c[d[(i, j)], d[(k, l)]] for i = 1:dim, j = 1:dim, k = 1:dim, l = 1:dim
    ])
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
for T in (:TensorStiffness, :TensorCompliance)
    eval(quote
        $T(a::AbstractArray) = $T(SArray{Tuple{3,3,3,3}}(a))
    end)
end

const VOIGT_INDICES =
    ((1, 1), (2, 2), (3, 3), (3, 2), (3, 1), (2, 1), (2, 3), (1, 3), (1, 2))

include("StabilityCriteria.jl")
include("Moduli.jl")

end
