module Moduli

export BulkModulus,
    YoungModulus, Lame1stParameter, ShearModulus, PoissonRatio, LongitudinalModulus

abstract type ElasticModulus{T} end
struct BulkModulus{T} <: ElasticModulus{T}
    v::T
end
struct YoungModulus{T} <: ElasticModulus{T}
    v::T
end
struct Lame1stParameter{T} <: ElasticModulus{T}
    v::T
end
struct ShearModulus{T} <: ElasticModulus{T}
    v::T
end
struct PoissonRatio{T} <: ElasticModulus{T}
    v::T
end
struct LongitudinalModulus{T} <: ElasticModulus{T}
    v::T
end

BulkModulus(E::YoungModulus, λ::Lame1stParameter) = BulkModulus((E.v + 3λ.v + _auxiliaryR(E, λ)) / 6)
BulkModulus(E::YoungModulus, G::ShearModulus) = BulkModulus(E.v * G.v / (3 * (3G.v - E.v)))
BulkModulus(E::YoungModulus, ν::PoissonRatio) = BulkModulus(E.v / (3 * (1 - 2ν.v)))
BulkModulus(E::YoungModulus, M::LongitudinalModulus) = BulkModulus((3M.v - E.v + _auxiliaryS(E, M)))
BulkModulus(λ::Lame1stParameter, G::ShearModulus) = BulkModulus(λ.v + 2G.v / 3)
BulkModulus(λ::Lame1stParameter, ν::PoissonRatio) = BulkModulus(λ.v * (1 + ν.v) / 3 / ν.v)
BulkModulus(λ::Lame1stParameter, M::LongitudinalModulus) = BulkModulus((M.v + 2λ.v) / 3)
BulkModulus(G::ShearModulus, ν::PoissonRatio) = BulkModulus(2G.v * (1 + ν.v) / (3 * (1 - 2ν.v)))
BulkModulus(G::ShearModulus, M::LongitudinalModulus) = BulkModulus(M.v - 4G.v / 3)
BulkModulus(ν::PoissonRatio, M::LongitudinalModulus) = BulkModulus(M.v * (1 + ν.v) / (3 * (1 - ν.v)))
BulkModulus(x::ElasticModulus, y::ElasticModulus) = BulkModulus(y, x)

YoungModulus(K::BulkModulus, λ::Lame1stParameter) = YoungModulus(9K.v * (K.v - λ.v) / (3K.v - λ.v))
YoungModulus(K::BulkModulus, G::ShearModulus) = YoungModulus(9K.v * G.v / (3K.v + G.v))
YoungModulus(K::BulkModulus, ν::PoissonRatio) = YoungModulus(3K.v * (1 - 2ν.v))
YoungModulus(K::BulkModulus, M::LongitudinalModulus) = YoungModulus(9K.v * (M.v - K.v) / (3K.v + M.v))
YoungModulus(λ::Lame1stParameter, G::ShearModulus) = YoungModulus(G.v * (3λ.v + 2G.v) / (λ.v + G.v))
YoungModulus(λ::Lame1stParameter, ν::PoissonRatio) = YoungModulus(λ.v * (1 + ν.v) * (1 - 2ν.v) / ν.v)
YoungModulus(λ::Lame1stParameter, M::LongitudinalModulus) = YoungModulus((M.v - λ.v) * (M.v + 2λ.v) / (M.v + λ.v))
YoungModulus(G::ShearModulus, ν::PoissonRatio) = YoungModulus(2G.v * (1 + ν.v))
YoungModulus(G::ShearModulus, M::LongitudinalModulus) = YoungModulus(G.v * (3M.v - 4G.v) / (M.v - G.v))
YoungModulus(ν::PoissonRatio, M::LongitudinalModulus) = YoungModulus(M.v * (1 + ν.v) * (1 - 2ν.v) / (1 - ν.v))

# These are helper functions and should not be exported!
_auxiliaryR(E::YoungModulus, λ::Lame1stParameter) = sqrt(E.v^2 + 9 * λ.v^2 + 2E.v * λ.v)
_auxiliaryS(E::YoungModulus, M::LongitudinalModulus) = sqrt(E.v^2 + 9 * M.v^2 - 10E.v * M.v)

for T in (
    :BulkModulus,
    :YoungModulus,
    :Lame1stParameter,
    :ShearModulus,
    :PoissonRatio,
    :LongitudinalModulus,
)
    eval(quote
        # As long as creating from one parameter already defined, return the parameter itself.
        $T(x::$T, ::ElasticModulus) = x
        # If there is no method match, switch the parameters and try again.
        $T(x::ElasticModulus, y::ElasticModulus) = $T(y, x)
    end)
end

end # module Moduli
