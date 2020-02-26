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

# These are helper functions and should not be exported!
_auxiliaryR(E::YoungModulus, λ::Lame1stParameter) = sqrt(E.v^2 + 9 * λ.v^2 + 2E.v * λ.v)
_auxiliaryS(E::YoungModulus, M::LongitudinalModulus) = sqrt(E.v^2 + 9 * M.v^2 - 10E.v * M.v)

BulkModulus(K::BulkModulus, ::ElasticModulus) = K
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

YoungModulus(E::YoungModulus, ::ElasticModulus) = E

end # module Moduli
