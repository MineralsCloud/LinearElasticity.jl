module Isotropic

export bulk, young, lame1st, shear, poisson, longitudinal

bulk(; kwargs...) = bulk(NamedTuple(kwargs))
bulk((E, λ)::NamedTuple{(:E, :λ)}) = (E + 3 * λ + _R(E, λ)) / 6
bulk((E, G)::NamedTuple{(:E, :G)}) = E * G / (3 * (3 * G - E))
bulk((E, ν)::NamedTuple{(:E, :ν)}) = E / (3 * (1 - 2 * ν))
bulk((E, M)::NamedTuple{(:E, :M)}) = (3 * M - E + _S(E, M)) / 6
bulk((λ, G)::NamedTuple{(:λ, :G)}) = λ + 2 * G / 3
bulk((λ, ν)::NamedTuple{(:λ, :ν)}) = λ * (1 + ν) / 3 / ν
bulk((λ, M)::NamedTuple{(:λ, :M)}) = (M + 2 * λ) / 3
bulk((G, ν)::NamedTuple{(:G, :ν)}) = 2 * G * (1 + ν) / (3 * (1 - 2 * ν))
bulk((G, M)::NamedTuple{(:G, :M)}) = M - 4 * G / 3
bulk((ν, M)::NamedTuple{(:ν, :M)}) = M * (1 + ν) / (3 * (1 - ν))
bulk(x::NamedTuple) = haskey(x, :K) ? x[:K] : bulk(_reverse(x))

young(; kwargs...) = young(NamedTuple(kwargs))
young((K, λ)::NamedTuple{(:K, :λ)}) = 9 * K * (K - λ) / (3 * K - λ)
young((K, G)::NamedTuple{(:K, :G)}) = 9 * K * G / (3 * K + G)
young((K, ν)::NamedTuple{(:K, :ν)}) = 3 * K * (1 - 2 * ν)
young((K, M)::NamedTuple{(:K, :M)}) = 9 * K * (M - K) / (3 * K + M)
young((λ, G)::NamedTuple{(:λ, :G)}) = G * (3 * λ + 2 * G) / (λ + G)
young((λ, ν)::NamedTuple{(:λ, :ν)}) = λ * (1 + ν) * (1 - 2 * ν) / ν
young((λ, M)::NamedTuple{(:λ, :M)}) = (M - λ) * (M + 2 * λ) / (M + λ)
young((G, ν)::NamedTuple{(:G, :ν)}) = 2 * G * (1 + ν)
young((G, M)::NamedTuple{(:G, :M)}) = G * (3 * M - 4 * G) / (M - G)
young((ν, M)::NamedTuple{(:ν, :M)}) = M * (1 + ν) * (1 - 2 * ν) / (1 - ν)
young(x::NamedTuple) = haskey(x, :E) ? x[:E] : young(_reverse(x))

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
lamé1st(x::NamedTuple) = haskey(x, :λ) ? x[:λ] : lamé1st(_reverse(x))
const lame1st = lamé1st

shear(; kwargs...) = shear(NamedTuple(kwargs))
shear((K, E)::NamedTuple{(:K, :E)}) = 3K * E / (9K - E)
shear((K, λ)::NamedTuple{(:K, :λ)}) = 3(K - λ) / 2
shear((K, ν)::NamedTuple{(:K, :ν)}) = 3K * (1 - 2ν) / 2 / (1 + ν)
shear((K, M)::NamedTuple{(:K, :M)}) = 3(M - K) / 4
shear((E, λ)::NamedTuple{(:E, :λ)}) = (E - 3λ + _R(E, λ)) / 4
shear((E, ν)::NamedTuple{(:E, :ν)}) = E / 2 / (1 + ν)
shear((E, M)::NamedTuple{(:E, :M)}) = (3M + E - _S(E, M)) / 8
shear((λ, ν)::NamedTuple{(:λ, :ν)}) = λ * (1 - 2ν) / 2 / ν
shear((λ, M)::NamedTuple{(:λ, :M)}) = (M - λ) / 2
shear((ν, M)::NamedTuple{(:ν, :M)}) = M * (1 - 2ν) / 2 / (1 - ν)
shear(x::NamedTuple) = haskey(x, :G) ? x[:G] : shear(_reverse(x))
const lamé2nd = shear
const lame2nd = shear

poisson(; kwargs...) = poisson(NamedTuple(kwargs))
poisson((K, E)::NamedTuple{(:K, :E)}) = (3K - E) / 6 / K
poisson((K, λ)::NamedTuple{(:K, :λ)}) = λ / (3K - λ)
poisson((K, G)::NamedTuple{(:K, :G)}) = (3K - 2G) / 2 / (3K + G)
poisson((K, M)::NamedTuple{(:K, :M)}) = (3K - M) / (3K + M)
poisson((E, λ)::NamedTuple{(:E, :λ)}) = 2λ / (E + λ + _R(E, λ))
poisson((E, G)::NamedTuple{(:E, :G)}) = E / 2 / G - 1
poisson((E, M)::NamedTuple{(:E, :M)}) = (E - M + _S(E, M)) / 4 / M
poisson((λ, G)::NamedTuple{(:λ, :G)}) = λ / 2 / (λ + G)
poisson((λ, M)::NamedTuple{(:λ, :M)}) = λ / (M + λ)
poisson((G, M)::NamedTuple{(:G, :M)}) = (M - 2G) / 2 / (M - G)
poisson(x::NamedTuple) = haskey(x, :ν) ? x[:ν] : poisson(_reverse(x))

longitudinal(; kwargs...) = longitudinal(NamedTuple(kwargs))
longitudinal((K, E)::NamedTuple{(:K, :E)}) = 3K * (3K + E) / (9K - E)
longitudinal((K, λ)::NamedTuple{(:K, :λ)}) = 3K - 2λ
longitudinal((K, G)::NamedTuple{(:K, :G)}) = K + 4G / 3
longitudinal((K, ν)::NamedTuple{(:K, :ν)}) = 3K * (1 - ν) / (1 + ν)
longitudinal((E, λ)::NamedTuple{(:E, :λ)}) = (E - λ + _R(E, λ)) / 2
longitudinal((E, G)::NamedTuple{(:E, :G)}) = G * (4G - E) / (3G - E)
longitudinal((E, ν)::NamedTuple{(:E, :ν)}) = E * (1 - ν) / (1 + ν) / (1 - 2ν)
longitudinal((λ, G)::NamedTuple{(:λ, :G)}) = λ + 2G
longitudinal((λ, ν)::NamedTuple{(:λ, :ν)}) = λ * (1 - ν) / ν
longitudinal((G, ν)::NamedTuple{(:G, :ν)}) = 2G * (1 - ν) / (1 - 2ν)
longitudinal(x::NamedTuple) = haskey(x, :M) ? x[:M] : longitudinal(_reverse(x))
const constrained = longitudinal

# These are helper functions and should not be exported!
_R(E, λ) = sqrt(E^2 + 9 * λ^2 + 2 * E * λ)
_S(E, M) = sqrt(E^2 + 9 * M^2 - 10 * E * M)

_reverse(x::NamedTuple) = (; zip(reverse(propertynames(x)), reverse(values(x)))...)

end
