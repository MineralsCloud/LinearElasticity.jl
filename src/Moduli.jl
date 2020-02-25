module Moduli

abstract type ElasticModulus{T} end
struct BulkModulus{T}
    v::T
end
struct YoungModulus{T}
    v::T
end
struct Lame1stParameter{T}
    v::T
end
struct ShearModulus{T}
    v::T
end
struct PoissonRatio{T}
    v::T
end
struct LongitudinalModulus{T}
    v::T
end

# These are helper functions and should not be exported!
_auxiliaryR(E::YoungModulus, λ::Lame1stParameter) = sqrt(E.v^2 + 9 * λ.v^2 + 2E.v * λ.v)

(::BulkModulus)(E::YoungModulus, λ::Lame1stParameter) = (E + 3λ + _auxiliaryR(E, λ)) / 6
(::BulkModulus)(E::YoungModulus, G::ShearModulus) = E * G / (3 * (3G - E))

end # module Moduli
