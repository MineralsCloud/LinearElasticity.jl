module Isotropic

export BulkModulus,
    YoungModulus,
    Lamé1stParameter,
    Lame1stParameter,
    ShearModulus,
    Lamé2ndParameter,
    Lame2ndParameter,
    PoissonRatio,
    LongitudinalModulus,
    unpack

abstract type ElasticModulus{T} end
struct BulkModulus{T} <: ElasticModulus{T}
    v::T
end
struct YoungModulus{T} <: ElasticModulus{T}
    v::T
end
struct Lamé1stParameter{T} <: ElasticModulus{T}
    v::T
end
const Lame1stParameter = Lamé1stParameter
struct ShearModulus{T} <: ElasticModulus{T}
    v::T
end
const Lamé2ndParameter = ShearModulus
const Lame2ndParameter = Lamé2ndParameter
struct PoissonRatio{T} <: ElasticModulus{T}
    v::T
end
struct LongitudinalModulus{T} <: ElasticModulus{T}
    v::T
end

BulkModulus(E::YoungModulus, λ::Lamé1stParameter) =
    BulkModulus((E.v + 3 * λ.v + _auxiliaryR(E, λ)) / 6)
BulkModulus(E::YoungModulus, G::ShearModulus) =
    BulkModulus(E.v * G.v / (3 * (3 * G.v - E.v)))
BulkModulus(E::YoungModulus, ν::PoissonRatio) = BulkModulus(E.v / (3 * (1 - 2 * ν.v)))
BulkModulus(E::YoungModulus, M::LongitudinalModulus) =
    BulkModulus((3 * M.v - E.v + _auxiliaryS(E, M)))
BulkModulus(λ::Lamé1stParameter, G::ShearModulus) = BulkModulus(λ.v + 2 * G.v / 3)
BulkModulus(λ::Lamé1stParameter, ν::PoissonRatio) = BulkModulus(λ.v * (1 + ν.v) / 3 / ν.v)
BulkModulus(λ::Lamé1stParameter, M::LongitudinalModulus) = BulkModulus((M.v + 2 * λ.v) / 3)
BulkModulus(G::ShearModulus, ν::PoissonRatio) =
    BulkModulus(2 * G.v * (1 + ν.v) / (3 * (1 - 2 * ν.v)))
BulkModulus(G::ShearModulus, M::LongitudinalModulus) = BulkModulus(M.v - 4 * G.v / 3)
BulkModulus(ν::PoissonRatio, M::LongitudinalModulus) =
    BulkModulus(M.v * (1 + ν.v) / (3 * (1 - ν.v)))

YoungModulus(K::BulkModulus, λ::Lamé1stParameter) =
    YoungModulus(9 * K.v * (K.v - λ.v) / (3 * K.v - λ.v))
YoungModulus(K::BulkModulus, G::ShearModulus) =
    YoungModulus(9 * K.v * G.v / (3 * K.v + G.v))
YoungModulus(K::BulkModulus, ν::PoissonRatio) = YoungModulus(3 * K.v * (1 - 2 * ν.v))
YoungModulus(K::BulkModulus, M::LongitudinalModulus) =
    YoungModulus(9 * K.v * (M.v - K.v) / (3 * K.v + M.v))
YoungModulus(λ::Lamé1stParameter, G::ShearModulus) =
    YoungModulus(G.v * (3 * λ.v + 2 * G.v) / (λ.v + G.v))
YoungModulus(λ::Lamé1stParameter, ν::PoissonRatio) =
    YoungModulus(λ.v * (1 + ν.v) * (1 - 2 * ν.v) / ν.v)
YoungModulus(λ::Lamé1stParameter, M::LongitudinalModulus) =
    YoungModulus((M.v - λ.v) * (M.v + 2 * λ.v) / (M.v + λ.v))
YoungModulus(G::ShearModulus, ν::PoissonRatio) = YoungModulus(2 * G.v * (1 + ν.v))
YoungModulus(G::ShearModulus, M::LongitudinalModulus) =
    YoungModulus(G.v * (3 * M.v - 4 * G.v) / (M.v - G.v))
YoungModulus(ν::PoissonRatio, M::LongitudinalModulus) =
    YoungModulus(M.v * (1 + ν.v) * (1 - 2 * ν.v) / (1 - ν.v))

Lamé1stParameter(K::BulkModulus, E::YoungModulus) =
    Lamé1stParameter((9K.v^2 - 3K.v * E.v) / (9K.v - E.v))
Lamé1stParameter(K::BulkModulus, G::ShearModulus) = Lamé1stParameter(K.v - 2G.v / 3)
Lamé1stParameter(K::BulkModulus, ν::PoissonRatio) = Lamé1stParameter(3K.v * ν.v / (1 + ν.v))
Lamé1stParameter(K::BulkModulus, M::LongitudinalModulus) =
    Lamé1stParameter((3K.v - M.v) / 2)
Lamé1stParameter(E::YoungModulus, G::ShearModulus) =
    Lamé1stParameter(G.v * (E.v - 2G.v) / (3G.v - E.v))
Lamé1stParameter(E::YoungModulus, ν::PoissonRatio) =
    Lamé1stParameter(E.v * ν.v / (1 + ν.v) / (1 - 2ν.v))
Lamé1stParameter(E::YoungModulus, M::LongitudinalModulus) =
    Lamé1stParameter((M.v - E.v + _auxiliaryS(E, M)) / 4)
Lamé1stParameter(G::ShearModulus, ν::PoissonRatio) =
    Lamé1stParameter(2G.v * ν.v / (1 - 2ν.v))
Lamé1stParameter(G::ShearModulus, M::ShearModulus) = Lamé1stParameter(M.v - 2G.v)
Lamé1stParameter(ν::PoissonRatio, M::LongitudinalModulus) =
    Lamé1stParameter(M.v * ν.v / (1 - ν.v))

# These are helper functions and should not be exported!
_auxiliaryR(E::YoungModulus, λ::Lamé1stParameter) = sqrt(E.v^2 + 9 * λ.v^2 + 2 * E.v * λ.v)
_auxiliaryS(E::YoungModulus, M::LongitudinalModulus) =
    sqrt(E.v^2 + 9 * M.v^2 - 10 * E.v * M.v)

unpack(x::ElasticModulus) = getfield(x, :v)
unpack(x::ElasticModulus...) = unpack.(x)
unpack(x::AbstractArray{<:ElasticModulus}) = unpack.(x)

for T in (
    :BulkModulus,
    :YoungModulus,
    :Lame1stParameter,
    :ShearModulus,
    :PoissonRatio,
    :LongitudinalModulus,
)
    @eval begin  # See https://github.com/JuliaLang/julia/blob/45c518d/base/fastmath.jl#L253-L259
        # As long as creating from one parameter already defined, return the parameter itself.
        $T(x::$T, ::ElasticModulus) = x
        # If there is no method match, switch the parameters and try again.
        $T(x::ElasticModulus, y::ElasticModulus) = $T(y, x)
    end
end

end
