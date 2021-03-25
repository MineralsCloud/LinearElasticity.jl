module StabilityConditions

using LinearAlgebra: diag, issymmetric

using Crystallography:
    CrystalSystem, Cubic, Hexagonal, Tetragonal, Trigonal, Orthorhombic, Monoclinic

using LinearElasticity: EngineeringStiffness, EngineeringCompliance

export isstable, issystem

ispositive(x) = x > zero(x)

function criteria(::Cubic, c::EngineeringStiffness)
    c₁₁, c₁₂, c₄₄ = c[1, 1], c[1, 2], c[4, 4]
    return (  # Must satisfy all criteria!
        c₁₁ > abs(c₁₂),
        c₁₁ > -2c₁₂,
        ispositive(c₄₄),
    )
end
function criteria(::Hexagonal, c::EngineeringStiffness)
    c₁₁, c₁₂, c₁₃, c₃₃, c₄₄, c₆₆ = c[1, 1], c[1, 2], c[1, 3], c[3, 3], c[4, 4], c[6, 6]
    return (  # Must satisfy all criteria!
        2c₆₆ == c₁₁ - c₁₂,
        c₁₁ > abs(c₁₂),
        2c₁₃^2 < c₃₃ * (c₁₁ + c₁₂),
        ispositive(c₄₄),
        ispositive(c₆₆),
    )
end
function criteria(::Tetragonal, c::EngineeringStiffness)
    c₁₁, c₁₂, c₁₃, c₁₆, c₃₃, c₄₄, c₆₆ =
        c[1, 1], c[1, 2], c[1, 3], c[1, 6], c[3, 3], c[4, 4], c[6, 6]
    if iszero(c₁₆)  # Tetragonal (I) class
        return (
            c₁₁ > abs(c₁₂),
            2c₁₃^2 < c₃₃ * (c₁₁ + c₁₂),
            ispositive(c₄₄),
            ispositive(c₆₆),
        )
    else
        return (  # Tetragonal (II) class
            c₁₁ > abs(c₁₂),
            2c₁₃^2 < c₃₃ * (c₁₁ + c₁₂),
            ispositive(c₄₄),
            2c₁₆^2 < c₆₆ * (c₁₁ - c₁₂),
        )
    end
end
function criteria(::Trigonal, c::EngineeringStiffness)
    c₁₁, c₁₂, c₁₃, c₁₄, c₁₅, c₃₃, c₄₄, c₆₆ =
        c[1, 1], c[1, 2], c[1, 3], c[1, 4], c[1, 5], c[3, 3], c[4, 4], c[6, 6]
    return (
        c₁₁ > abs(c₁₂),
        ispositive(c₄₄),
        2c₁₃^2 < c₃₃ * (c₁₁ + c₁₂),
        c₄₄ * (c₁₁ - c₁₂) == 2c₄₄ * c₆₆,
        if iszero(c₁₅)  # Rhombohedral (I) class
            2c₁₄^2 < c₄₄ * (c₁₁ - c₁₂)
        else  # Rhombohedral (II) class
            2(c₁₄^2 + c₁₅^2) < c₄₄ * (c₁₁ - c₁₂)
        end,
    )
end
function criteria(::Orthorhombic, c::EngineeringStiffness)
    c₁₁, c₂₂, c₃₃, c₄₄, c₅₅, c₆₆ = diag(c)
    c₁₂, c₁₃, c₂₃ = c[1, 2], c[1, 3], c[2, 3]
    return (
        map(ispositive, (c₁₁, c₄₄, c₅₅, c₆₆))...,
        c₁₁ * c₂₂ > c₁₂^2,
        c₁₁ * c₂₂ * c₃₃ + 2c₁₂ * c₁₃ * c₂₃ > c₁₁ * c₂₃^2 + c₂₂ * c₁₃^2 + c₃₃ * c₁₂^2,
    )
end
function criteria(::Monoclinic, c::EngineeringStiffness)
    c₁₁, c₂₂, c₃₃, c₄₄, c₅₅, c₆₆ = diag(c)
    c₁₂, c₁₃, c₁₅, c₂₃, c₂₅, c₃₅, c₄₆ =
        c[1, 2], c[1, 3], c[1, 5], c[2, 3], c[2, 5], c[3, 5], c[4, 6]
    g =
        c₁₁ * c₂₂ * c₃₃ - c₁₁ * c₂₃ * c₂₃ - c₂₂ * c₁₃ * c₁₃ - c₃₃ * c₁₂ * c₁₂ +
        2c₁₂ * c₁₃ * c₂₃
    return (
        map(ispositive, (c₁₁, c₂₂, c₃₃, c₄₄, c₅₅, c₆₆))...,
        c₁₁ + c₂₂ + c₃₃ > -2(c₁₂ + c₁₃ + c₂₃),
        c₃₃ * c₅₅ > c₃₅^2,
        c₄₄ * c₆₆ > c₄₆^2,
        c₂₂ + c₃₃ > 2c₂₃,
        c₂₂ * (c₃₃ * c₅₅ - c₃₅^2) + 2c₂₃ * c₂₅ * c₃₅ - c₂₃^2 * c₅₅ > c₂₅^2 * c₃₃,
        2(
            c₁₅ * c₂₅ * (c₃₃ * c₁₂ - c₁₃ * c₂₃) +
            c₁₅ * c₃₅ * (c₂₂ * c₁₃ - c₁₂ * c₂₃) +
            c₂₅ * c₃₅ * (c₁₁ * c₂₃ - c₁₂ * c₁₃)
        ) - (
            c₁₅ * c₁₅ * (c₂₂ * c₃₃ - c₂₃^2) +
            c₂₅ * c₂₅ * (c₁₁ * c₃₃ - c₁₃^2) +
            c₃₅ * c₃₅ * (c₁₁ * c₂₂ - c₁₂^2)
        ) > -c₅₅ * g,
    )
end
criteria(C::CrystalSystem, s::EngineeringCompliance) = criteria(C, inv(s))

isstable(C::CrystalSystem, x::Union{EngineeringStiffness,EngineeringCompliance}) =
    all(criteria(C, x))

function issystem(::Cubic, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        c[1, 1] == c[2, 2] == c[3, 3],
        c[4, 4] == c[5, 5] == c[6, 6],
        c[1, 2] == c[1, 3] == c[2, 3],
    ))
end
function issystem(::Hexagonal, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        2c[6, 6] == c[1, 1] - c[1, 2],
    ))
end
function issystem(::Tetragonal, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        if iszero(c[1, 6])  # Tetragonal (I) class, 4mm, -42m, 422, 4/mmm
            true
        else  # Tetragonal (II) class, 4, -4, 4/m
            c[1, 6] == -c[2, 6]
        end,
    ))
end
function issystem(::Trigonal, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        2c[6, 6] == c[1, 1] - c[1, 2],
        c[1, 4] == -c[2, 4] == -c[5, 6],
        if iszero(c[1, 5])  # # Rhombohedral (I) class, 32, -3m, 3m
            true
        else  # Rhombohedral (II) class, 3, -3
            -c[1, 5] == c[2, 5] == c[4, 6]
        end,
    ))
end
function issystem(::Orthorhombic, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        all(
            iszero,
            (
                c[1, 4],
                c[2, 4],
                c[3, 4],
                c[1, 5],
                c[2, 5],
                c[3, 5],
                c[4, 5],
                c[1, 6],
                c[2, 6],
                c[3, 6],
                c[4, 6],
                c[5, 6],
            ),
        ),
        all(
            !iszero,
            (
                c[1, 1],
                c[1, 2],
                c[1, 3],
                c[2, 2],
                c[2, 3],
                c[3, 3],
                c[4, 4],
                c[5, 5],
                c[6, 6],
            ),
        ),
    ))
end
function issystem(::Monoclinic, c::EngineeringStiffness)
    return all((
        issymmetric(c),
        all(
            iszero,
            (c[1, 4], c[2, 4], c[3, 4], c[4, 5], c[1, 6], c[2, 6], c[3, 6], c[5, 6]),
        ),
    ))
end
issystem(C::CrystalSystem, s::EngineeringCompliance) = issystem(C, inv(s))

end
