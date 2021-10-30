Base.inv(c::TensorStiffness) = TensorCompliance(inv(c))
Base.inv(s::TensorCompliance) = TensorStiffness(inv(s))
Base.inv(c::EngineeringStiffness) = EngineeringCompliance(inv(c))
Base.inv(s::EngineeringCompliance) = EngineeringStiffness(inv(s))

const VOIGT_INDICES =
    ((1, 1), (2, 2), (3, 3), (3, 2), (3, 1), (2, 1), (2, 3), (1, 3), (1, 2))

Base.convert(::Type{TensorStress{T}}, s::EngineeringStress{T}) where {T} =
    TensorStress([s[1], s[6], s[5], s[2], s[4], s[3]])
Base.convert(::Type{EngineeringStress{T}}, σ::TensorStress{T}) where {T} =
    EngineeringStress([σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2]])
Base.convert(::Type{TensorStrain{T}}, ϵ::EngineeringStrain{T}) where {T} =
    TensorStrain([ϵ[1], ϵ[6] / 2, ϵ[5] / 2, ϵ[2], ϵ[4] / 2, ϵ[3]])
Base.convert(::Type{EngineeringStrain{T}}, ε::TensorStrain{T}) where {T} =
    EngineeringStrain([ε[1, 1], ε[2, 2], ε[3, 3], 2ε[2, 3], 2ε[1, 3], 2ε[1, 2]])
function Base.convert(::Type{EngineeringStiffness{T}}, c::TensorStiffness{T}) where {T}
    p, dim = pairs(VOIGT_INDICES), 6
    # From https://github.com/KristofferC/Tensors.jl/blob/bff451c/src/utilities.jl#L5-L14
    return EngineeringStiffness([c[p[i]..., p[j]...] for i in 1:dim, j in 1:dim])
end
function Base.convert(::Type{TensorStiffness{T}}, c::EngineeringStiffness{T}) where {T}
    d, dim = Dict(zip(VOIGT_INDICES, [1, 2, 3, 4, 5, 6, 4, 5, 6])), 3
    return TensorStiffness([
        c[d[(i, j)], d[(k, l)]] for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim
    ])
end
function Base.convert(::Type{EngineeringCompliance{T}}, s::TensorCompliance{T}) where {T}  # FIXME Rules are wrong
    p, dim = pairs(VOIGT_INDICES), 6
    # From https://github.com/KristofferC/Tensors.jl/blob/bff451c/src/utilities.jl#L5-L14
    return EngineeringStiffness([s[p[i]..., p[j]...] for i in 1:dim, j in 1:dim])
end
function Base.convert(::Type{TensorCompliance{T}}, s::EngineeringCompliance{T}) where {T}  # FIXME Rules are wrong
    d, dim = Dict(zip(VOIGT_INDICES, [1, 2, 3, 4, 5, 6, 4, 5, 6])), 3
    return TensorStiffness([
        s[d[(i, j)], d[(k, l)]] for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim
    ])
end
