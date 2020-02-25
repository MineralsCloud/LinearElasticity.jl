module Moduli

export BulkModulus,
    YoungModulus,
    Lame1stParameter,
    ShearModulus,
    PoissonRatio,
    LongitudinalModulus

abstract type ElasticModulus{T} end
struct BulkModulus{T} <: ElasticModulus{T} end
struct YoungModulus{T} <: ElasticModulus{T} end
struct Lame1stParameter{T} <: ElasticModulus{T} end
struct ShearModulus{T} <: ElasticModulus{T} end
struct PoissonRatio{T} <: ElasticModulus{T} end
struct LongitudinalModulus{T} <: ElasticModulus{T} end

# These are helper functions and should not be exported!
_auxiliaryR(E, λ) = sqrt(E^2 + 9 * λ^2 + 2E * λ)
_auxiliaryS(E, M) = sqrt(E^2 + 9 * M^2 - 10E * M)

BulkModulus(::YoungModulus{E}, ::Lame1stParameter{λ}) where {E,λ} =
    BulkModulus((E + 3λ + _auxiliaryR(E, λ)) / 6)
BulkModulus(::YoungModulus{E}, ::ShearModulus{G}) where {E,G} =
    BulkModulus(E * G / (3 * (3G - E)))

for T in (
    :BulkModulus,
    :YoungModulus,
    :Lame1stParameter,
    :ShearModulus,
    :PoissonRatio,
    :LongitudinalModulus,
)
    eval(quote
        $T(v) = $T{v}()
        (::$T{v})() where {v} = v
    end)
end

end # module Moduli
