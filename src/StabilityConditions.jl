module StabilityConditions

using LinearAlgebra: diag

using Crystallography:
    CrystalSystem, Cubic, Hexagonal, Tetragonal, Trigonal, Orthorhombic, Monoclinic

using LinearElasticity: EngineeringStiffness, EngineeringCompliance

export isstable

function isstable(::Cubic, c::EngineeringStiffness)
    c11, c12, c44 = c[1, 1], c[1, 2], c[1, 4]
    return all([  # Must satisfy all criteria!
        c11 > abs(c12),
        c11 + 2 * c12 > 0,
        c44 > 0,
    ])
end # function isstable
function isstable(::Hexagonal, c::EngineeringStiffness)
    c11, c12, c13, c33, c44, c66 = c[1, 1], c[1, 2], c[1, 3], c[3, 3], c[4, 4], c[6, 6]
    return all([  # Must satisfy all criteria!
        c66 == (c11 - c12) / 2,
        c11 > abs(c12),
        2 * c13^2 < c33 * (c11 + c12),
        c44 > 0,
        c66 > 0,
    ])
end # function isstable
function isstable(::Tetragonal, c::EngineeringStiffness)
    c11, c12, c13, c16, c33, c44, c66 =
        c[1, 1], c[1, 2], c[1, 3], c[1, 6], c[3, 3], c[4, 4], c[6, 6]
    if c16 == 0  # Tetragonal (I) class
        return all([c11 > abs(c12), 2 * c13^2 < c33 * (c11 + c12), c44 > 0, c66 > 0])
    end
    return all([  # Tetragonal (II) class
        c11 > abs(c12),
        2 * c13^2 < c33 * (c11 + c12),
        c44 > 0,
        2 * c16^2 < c66 * (c11 - c12),
    ])
end # function isstable
function isstable(::Trigonal, c::EngineeringStiffness)
    c11, c12, c13, c14, c15, c33, c44, c66 =
        c[1, 1], c[1, 2], c[1, 3], c[1, 4], c[1, 5], c[3, 3], c[4, 4], c[6, 6]
    if c15 == 0  # Rhombohedral (I) class
        return all([
            c11 > abs(c12),
            c44 > 0,
            c13^2 < 0.5 * c33 * (c11 + c12),
            c14^2 < 0.5 * c44 * (c11 - c12),
            0.5 * c44 * (c11 - c12) == c44 * c66,
        ])
    end
    return all([  # Rhombohedral (II) class
        c11 > abs(c12),
        c44 > 0,
        c13^2 < 0.5 * c33 * (c11 + c12),
        c14^2 + c15^2 < 0.5 * c44 * (c11 - c12),
        0.5 * c44 * (c11 - c12) == c44 * c66,
    ])
end # function isstable
function isstable(::Orthorhombic, c::EngineeringStiffness)
    c11, c22, c33, c44, c55, c66 = diag(c)
    c12, c13, c23 = c[1, 2], c[1, 3], c[2, 3]
    return all([
        c11 > 0,
        c11 * c22 > c12^2,
        c11 * c22 * c33 + 2 * c12 * c13 * c23 > c11 * c23^2 + c22 * c13^2 + c33 * c12^2,
        c44 > 0,
        c55 > 0,
        c66 > 0,
    ])
end # function isstable
function isstable(::Monoclinic, c::EngineeringStiffness)
    c11, c22, c33, c44, c55, c66 = diag(c)
    c12, c13, c15, c23, c25, c35, c46 =
        c[1, 2], c[1, 3], c[1, 5], c[2, 3], c[2, 5], c[3, 5], c[4, 6]
    g = c11 * c22 * c33 - c11 * c23 * c23 - c22 * c13 * c13 - c33 * c12 * c12 +
        2 * c12 * c13 * c23
    return all([
        all([c11, c22, c33, c44, c55, c66] .> 0),
        c11 + c22 + c33 + 2 * (c12 + c13 + c23) > 0,
        c33 * c55 - c35^2 > 0,
        c44 * c66 - c46^2 > 0,
        c22 + c33 - 2 * c23 > 0,
        c22 * (c33 * c55 - c35^2) + 2 * c23 * c25 * c35 - c23^2 * c55 - c25^2 * c33 > 0,
        2 * (
            c15 * c25 * (c33 * c12 - c13 * c23) +
            c15 * c35 * (c22 * c13 - c12 * c23) +
            c25 * c35 * (c11 * c23 - c12 * c13)
        ) - (
            c15 * c15 * (c22 * c33 - c23^2) +
            c25 * c25 * (c11 * c33 - c13^2) +
            c35 * c35 * (c11 * c22 - c12^2)
        ) + c55 * g > 0,
    ])
end # function isstable
isstable(C::CrystalSystem, s::EngineeringCompliance) = isstable(C, inv(s))

end # module StabilityConditions
