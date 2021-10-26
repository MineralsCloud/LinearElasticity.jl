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
    return EngineeringStiffness([c[p[i]..., p[j]...] for i in 1:dim, j in 1:dim])
end # function Base.convert
function Base.convert(::Type{TensorStiffness{T}}, c::EngineeringStiffness{T}) where {T}
    d, dim = Dict(zip(VOIGT_INDICES, [1, 2, 3, 4, 5, 6, 4, 5, 6])), 3
    return TensorStiffness([
        c[d[(i, j)], d[(k, l)]] for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim
    ])
end # function Base.convert
function Base.convert(::Type{EngineeringCompliance{T}}, c::TensorCompliance{T}) where {T}  # FIXME Rules are wrong
    p, dim = pairs(VOIGT_INDICES), 6
    # From https://github.com/KristofferC/Tensors.jl/blob/bff451c/src/utilities.jl#L5-L14
    return EngineeringStiffness([c[p[i]..., p[j]...] for i in 1:dim, j in 1:dim])
end # function Base.convert
function Base.convert(::Type{TensorCompliance{T}}, c::EngineeringCompliance{T}) where {T}  # FIXME Rules are wrong
    d, dim = Dict(zip(VOIGT_INDICES, [1, 2, 3, 4, 5, 6, 4, 5, 6])), 3
    return TensorStiffness([
        c[d[(i, j)], d[(k, l)]] for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim
    ])
end # function Base.convert
