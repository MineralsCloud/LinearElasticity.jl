module StabilityConditions

using LinearAlgebra: diag

using Crystallography:
    CrystalSystem, Cubic, Hexagonal, Tetragonal, Trigonal, Orthorhombic, Monoclinic

using LinearElasticity: EngineeringStiffness, EngineeringCompliance

export isstable

function isstable(::Cubic, c::EngineeringStiffness)
    c₁₁, c₁₂, c₄₄ = c[1, 1], c[1, 2], c[4, 4]
    return all((  # Must satisfy all criteria!
        c₁₁ > abs(c₁₂),
        c₁₁ + 2c₁₂ > 0,
        c₄₄ > 0,
    ))
end
function isstable(::Hexagonal, c::EngineeringStiffness)
    c₁₁, c₁₂, c₁₃, c₃₃, c₄₄, c₆₆ = c[1, 1], c[1, 2], c[1, 3], c[3, 3], c[4, 4], c[6, 6]
    return all((  # Must satisfy all criteria!
        2c₆₆ == c₁₁ - c₁₂,
        c₁₁ > abs(c₁₂),
        2c₁₃^2 < c₃₃ * (c₁₁ + c₁₂),
        c₄₄ > 0,
        c₆₆ > 0,
    ))
end
function isstable(::Tetragonal, c::EngineeringStiffness)
    c₁₁, c₁₂, c₁₃, c₁₆, c₃₃, c₄₄, c₆₆ =
        c[1, 1], c[1, 2], c[1, 3], c[1, 6], c[3, 3], c[4, 4], c[6, 6]
    if c₁₆ == 0  # Tetragonal (I) class
        return all((c₁₁ > abs(c₁₂), 2 * c₁₃^2 < c₃₃ * (c₁₁ + c₁₂), c₄₄ > 0, c₆₆ > 0))
    end
    return all((  # Tetragonal (II) class
        c₁₁ > abs(c₁₂),
        2c₁₃^2 < c₃₃ * (c₁₁ + c₁₂),
        c₄₄ > 0,
        2c₁₆^2 < c₆₆ * (c₁₁ - c₁₂),
    ))
end
function isstable(::Trigonal, c::EngineeringStiffness)
    c₁₁, c₁₂, c₁₃, c₁₄, c₁₅, c₃₃, c₄₄, c₆₆ =
        c[1, 1], c[1, 2], c[1, 3], c[1, 4], c[1, 5], c[3, 3], c[4, 4], c[6, 6]
    if c₁₅ == 0  # Rhombohedral (I) class
        return all((
            c₁₁ > abs(c₁₂),
            c₄₄ > 0,
            2c₁₃^2 < c₃₃ * (c₁₁ + c₁₂),
            2c₁₄^2 < c₄₄ * (c₁₁ - c₁₂),
            c₄₄ * (c₁₁ - c₁₂) == 2c₄₄ * c₆₆,
        ))
    end
    return all((  # Rhombohedral (II) class
        c₁₁ > abs(c₁₂),
        c₄₄ > 0,
        2c₁₃^2 < c₃₃ * (c₁₁ + c₁₂),
        2(c₁₄^2 + c₁₅^2) < c₄₄ * (c₁₁ - c₁₂),
        c₄₄ * (c₁₁ - c₁₂) == 2c₄₄ * c₆₆,
    ))
end
function isstable(::Orthorhombic, c::EngineeringStiffness)
    c₁₁, c₂₂, c₃₃, c₄₄, c₅₅, c₆₆ = diag(c)
    c₁₂, c₁₃, c₂₃ = c[1, 2], c[1, 3], c[2, 3]
    return all((
        c₁₁ > 0,
        c₁₁ * c₂₂ > c₁₂^2,
        c₁₁ * c₂₂ * c₃₃ + 2c₁₂ * c₁₃ * c₂₃ > c₁₁ * c₂₃^2 + c₂₂ * c₁₃^2 + c₃₃ * c₁₂^2,
        c₄₄ > 0,
        c₅₅ > 0,
        c₆₆ > 0,
    ))
end
function isstable(::Monoclinic, c::EngineeringStiffness)
    c₁₁, c₂₂, c₃₃, c₄₄, c₅₅, c₆₆ = diag(c)
    c₁₂, c₁₃, c₁₅, c₂₃, c₂₅, c₃₅, c₄₆ =
        c[1, 2], c[1, 3], c[1, 5], c[2, 3], c[2, 5], c[3, 5], c[4, 6]
    g =
        c₁₁ * c₂₂ * c₃₃ - c₁₁ * c₂₃ * c₂₃ - c₂₂ * c₁₃ * c₁₃ - c₃₃ * c₁₂ * c₁₂ +
        2c₁₂ * c₁₃ * c₂₃
    return all((
        all((c₁₁, c₂₂, c₃₃, c₄₄, c₅₅, c₆₆) .> 0),
        c₁₁ + c₂₂ + c₃₃ + 2(c₁₂ + c₁₃ + c₂₃) > 0,
        c₃₃ * c₅₅ - c₃₅^2 > 0,
        c₄₄ * c₆₆ - c₄₆^2 > 0,
        c₂₂ + c₃₃ - 2c₂₃ > 0,
        c₂₂ * (c₃₃ * c₅₅ - c₃₅^2) + 2c₂₃ * c₂₅ * c₃₅ - c₂₃^2 * c₅₅ - c₂₅^2 * c₃₃ > 0,
        2(
            c₁₅ * c₂₅ * (c₃₃ * c₁₂ - c₁₃ * c₂₃) +
            c₁₅ * c₃₅ * (c₂₂ * c₁₃ - c₁₂ * c₂₃) +
            c₂₅ * c₃₅ * (c₁₁ * c₂₃ - c₁₂ * c₁₃)
        ) - (
            c₁₅ * c₁₅ * (c₂₂ * c₃₃ - c₂₃^2) +
            c₂₅ * c₂₅ * (c₁₁ * c₃₃ - c₁₃^2) +
            c₃₅ * c₃₅ * (c₁₁ * c₂₂ - c₁₂^2)
        ) + c₅₅ * g > 0,
    ))
end
isstable(C::CrystalSystem, s::EngineeringCompliance) = isstable(C, inv(s))

end
