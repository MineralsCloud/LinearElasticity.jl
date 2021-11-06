module Isotropic

export bulk, young, lame1st, shear, poisson, longitudinal

const ALLOWED_KEYS = (:K, :E, :λ, :G, :ν, :M)

bulk_modulus(; kwargs...) = bulk_modulus(NamedTuple(kwargs))
bulk_modulus((E, λ)::NamedTuple{(:E, :λ)}) = (E + 3λ + _R(E, λ)) / 6
bulk_modulus((E, G)::NamedTuple{(:E, :G)}) = E * G / 3(3G - E)
bulk_modulus((E, ν)::NamedTuple{(:E, :ν)}) = E / 3(1 - 2ν)
bulk_modulus((E, M)::NamedTuple{(:E, :M)}) = (3M - E + _S(E, M)) / 6
bulk_modulus((λ, G)::NamedTuple{(:λ, :G)}) = λ + 2G / 3
bulk_modulus((λ, ν)::NamedTuple{(:λ, :ν)}) = λ * (1 + ν) / 3ν
bulk_modulus((λ, M)::NamedTuple{(:λ, :M)}) = (M + 2λ) / 3
bulk_modulus((G, ν)::NamedTuple{(:G, :ν)}) = 2G * (1 + ν) / 3(1 - 2ν)
bulk_modulus((G, M)::NamedTuple{(:G, :M)}) = M - 4G / 3
bulk_modulus((ν, M)::NamedTuple{(:ν, :M)}) = M * (1 + ν) / 3(1 - ν)
function bulk_modulus(x::NamedTuple)
    _checkkeys(x)
    return haskey(x, :K) ? x[:K] : bulk_modulus(_reverse(x))
end
const bulk = bulk_modulus

young_modulus(; kwargs...) = young_modulus(NamedTuple(kwargs))
young_modulus((K, λ)::NamedTuple{(:K, :λ)}) = 9K * (K - λ) / (3K - λ)
young_modulus((K, G)::NamedTuple{(:K, :G)}) = 9K * G / (3K + G)
young_modulus((K, ν)::NamedTuple{(:K, :ν)}) = 3K * (1 - 2ν)
young_modulus((K, M)::NamedTuple{(:K, :M)}) = 9K * (M - K) / (3K + M)
young_modulus((λ, G)::NamedTuple{(:λ, :G)}) = G * (3λ + 2G) / (λ + G)
young_modulus((λ, ν)::NamedTuple{(:λ, :ν)}) = λ * (1 + ν) * (1 - 2ν) / ν
young_modulus((λ, M)::NamedTuple{(:λ, :M)}) = (M - λ) * (M + 2λ) / (M + λ)
young_modulus((G, ν)::NamedTuple{(:G, :ν)}) = 2G * (1 + ν)
young_modulus((G, M)::NamedTuple{(:G, :M)}) = G * (3M - 4G) / (M - G)
young_modulus((ν, M)::NamedTuple{(:ν, :M)}) = M * (1 + ν) * (1 - 2ν) / (1 - ν)
function young_modulus(x::NamedTuple)
    _checkkeys(x)
    return haskey(x, :E) ? x[:E] : young_modulus(_reverse(x))
end
const young = young_modulus

lamé1st(; kwargs...) = lamé1st(NamedTuple(kwargs))
lamé1st((K, E)::NamedTuple{(:K, :E)}) = (9K^2 - 3K * E) / (9K - E)
lamé1st((K, G)::NamedTuple{(:K, :G)}) = K - 2G / 3
lamé1st((K, ν)::NamedTuple{(:K, :ν)}) = 3K * ν / (1 + ν)
lamé1st((K, M)::NamedTuple{(:K, :M)}) = (3K - M) / 2
lamé1st((E, G)::NamedTuple{(:E, :G)}) = G * (E - 2G) / (3G - E)
lamé1st((E, ν)::NamedTuple{(:E, :ν)}) = E * ν / (1 + ν) / (1 - 2ν)
lamé1st((E, M)::NamedTuple{(:E, :M)}) = (M - E + _S(E, M)) / 4
lamé1st((G, ν)::NamedTuple{(:G, :ν)}) = 2G * ν / (1 - 2ν)
lamé1st((G, M)::NamedTuple{(:G, :M)}) = M - 2G
lamé1st((ν, M)::NamedTuple{(:ν, :M)}) = M * ν / (1 - ν)
function lamé1st(x::NamedTuple)
    _checkkeys(x)
    return haskey(x, :λ) ? x[:λ] : lamé1st(_reverse(x))
end
const lame1st = lamé1st

shear_modulus(; kwargs...) = shear_modulus(NamedTuple(kwargs))
shear_modulus((K, E)::NamedTuple{(:K, :E)}) = 3K * E / (9K - E)
shear_modulus((K, λ)::NamedTuple{(:K, :λ)}) = 3(K - λ) / 2
shear_modulus((K, ν)::NamedTuple{(:K, :ν)}) = 3K * (1 - 2ν) / 2(1 + ν)
shear_modulus((K, M)::NamedTuple{(:K, :M)}) = 3(M - K) / 4
shear_modulus((E, λ)::NamedTuple{(:E, :λ)}) = (E - 3λ + _R(E, λ)) / 4
shear_modulus((E, ν)::NamedTuple{(:E, :ν)}) = E / 2(1 + ν)
shear_modulus((E, M)::NamedTuple{(:E, :M)}) = (3M + E - _S(E, M)) / 8
shear_modulus((λ, ν)::NamedTuple{(:λ, :ν)}) = λ * (1 - 2ν) / 2ν
shear_modulus((λ, M)::NamedTuple{(:λ, :M)}) = (M - λ) / 2
shear_modulus((ν, M)::NamedTuple{(:ν, :M)}) = M * (1 - 2ν) / 2(1 - ν)
function shear_modulus(x::NamedTuple)
    _checkkeys(x)
    return haskey(x, :G) ? x[:G] : shear_modulus(_reverse(x))
end
const lamé2nd = shear_modulus
const lame2nd = shear_modulus
const shear = shear_modulus

poisson_ratio(; kwargs...) = poisson_ratio(NamedTuple(kwargs))
poisson_ratio((K, E)::NamedTuple{(:K, :E)}) = (3K - E) / 6K
poisson_ratio((K, λ)::NamedTuple{(:K, :λ)}) = λ / (3K - λ)
poisson_ratio((K, G)::NamedTuple{(:K, :G)}) = (3K - 2G) / 2(3K + G)
poisson_ratio((K, M)::NamedTuple{(:K, :M)}) = (3K - M) / (3K + M)
poisson_ratio((E, λ)::NamedTuple{(:E, :λ)}) = 2λ / (E + λ + _R(E, λ))
poisson_ratio((E, G)::NamedTuple{(:E, :G)}) = E / 2G - 1
poisson_ratio((E, M)::NamedTuple{(:E, :M)}) = (E - M + _S(E, M)) / 4M
poisson_ratio((λ, G)::NamedTuple{(:λ, :G)}) = λ / 2(λ + G)
poisson_ratio((λ, M)::NamedTuple{(:λ, :M)}) = λ / (M + λ)
poisson_ratio((G, M)::NamedTuple{(:G, :M)}) = (M - 2G) / 2(M - G)
function poisson_ratio(x::NamedTuple)
    _checkkeys(x)
    return haskey(x, :ν) ? x[:ν] : poisson_ratio(_reverse(x))
end
const poisson = poisson_ratio

longitudinal_modulus(; kwargs...) = longitudinal_modulus(NamedTuple(kwargs))
longitudinal_modulus((K, E)::NamedTuple{(:K, :E)}) = 3K * (3K + E) / (9K - E)
longitudinal_modulus((K, λ)::NamedTuple{(:K, :λ)}) = 3K - 2λ
longitudinal_modulus((K, G)::NamedTuple{(:K, :G)}) = K + 4G / 3
longitudinal_modulus((K, ν)::NamedTuple{(:K, :ν)}) = 3K * (1 - ν) / (1 + ν)
longitudinal_modulus((E, λ)::NamedTuple{(:E, :λ)}) = (E - λ + _R(E, λ)) / 2
longitudinal_modulus((E, G)::NamedTuple{(:E, :G)}) = G * (4G - E) / (3G - E)
longitudinal_modulus((E, ν)::NamedTuple{(:E, :ν)}) = E * (1 - ν) / (1 + ν) / (1 - 2ν)
longitudinal_modulus((λ, G)::NamedTuple{(:λ, :G)}) = λ + 2G
longitudinal_modulus((λ, ν)::NamedTuple{(:λ, :ν)}) = λ * (1 - ν) / ν
longitudinal_modulus((G, ν)::NamedTuple{(:G, :ν)}) = 2G * (1 - ν) / (1 - 2ν)
function longitudinal_modulus(x::NamedTuple)
    _checkkeys(x)
    return haskey(x, :M) ? x[:M] : longitudinal_modulus(_reverse(x))
end
const longitudinal = longitudinal_modulus

# These are helper functions and should not be exported!
_R(E, λ) = sqrt(E^2 + 9λ^2 + 2E * λ)
_S(E, M) = sqrt(E^2 + 9M^2 - 10E * M)  # FIXME: ±S

_reverse(x::NamedTuple) = (; zip(reverse(propertynames(x)), reverse(values(x)))...)

_checkkeys(x) = @assert all(key ∈ ALLOWED_KEYS for key in propertynames(x))

end
