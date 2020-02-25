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
_auxiliaryR(E::YoungModulus, λ::Lame1stParameter) = sqrt(E^2 + 9 * λ^2 + 2E * λ)
_auxiliaryS(E::YoungModulus, M::LongitudinalModulus) = sqrt(E^2 + 9 * M^2 - 10E * M)

BulkModulus(E::YoungModulus, λ::Lame1stParameter{}) =
    BulkModulus((E + 3λ + _auxiliaryR(E, λ)) / 6)
BulkModulus(E::YoungModulus, G::ShearModulus) = BulkModulus(E * G / (3 * (3G - E)))

for op in (:+, :-, :*, :/, ://, :div)
    eval(quote
        Base.$op(a::ElasticModulus, b::ElasticModulus) = $op(a.v, b.v)
        Base.$op(a::T, b::Number) where {T<:ElasticModulus} = T($op(a.v, b))
        Base.$op(a::Number, b::ElasticModulus) = $op(b, a)
    end)
end
for op in (:-,)
    eval(:(Base.$op(a::ElasticModulus) = $op(a.v)))
end
Base.:^(a::T, n::Number) where {T<:ElasticModulus} = T(a.v^n)

end # module Moduli
