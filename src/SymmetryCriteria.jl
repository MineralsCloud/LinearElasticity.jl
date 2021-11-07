module SymmetryCriteria

using CrystallographyBase:
    CrystalSystem,
    Cubic,
    Hexagonal,
    Tetragonal,
    Trigonal,
    Orthorhombic,
    Monoclinic,
    Triclinic
using ..LinearElasticity: StiffnessMatrix, ComplianceMatrix

export issystem, whichsystem

function symmetry_criteria(::Cubic, c::StiffnessMatrix)
    return (
        all(iszero, (c[1:3, 4:6]..., c[4, 5:6]..., c[5, 6])),
        all(!iszero, (c[1, 1], c[4, 4], c[1, 3])),
        c[1, 1] == c[2, 2] == c[3, 3],
        c[4, 4] == c[5, 5] == c[6, 6],
        c[1, 2] == c[1, 3] == c[2, 3],
    )
end
function symmetry_criteria(::Hexagonal, c::StiffnessMatrix)
    return (
        all(iszero, (c[1:3, 4:6]..., c[4, 5:6]..., c[5, 6])),
        all(!iszero, (c[1, 1], c[3, 3], c[4, 4], c[1, 2], c[1, 3], c[6, 6])),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        2c[6, 6] == c[1, 1] - c[1, 2],
    )
end
function symmetry_criteria(::Tetragonal, c::StiffnessMatrix)
    return (
        all(iszero, (c[1:3, 4:5]..., c[3, 6], c[4, 5:6]..., c[5, 6])),
        all(!iszero, (c[1, 1], c[1, 2], c[1, 3], c[3, 3], c[4, 4], c[6, 6])),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        if all(iszero, (c[2, 6], c[1, 6]))  # Tetragonal (I) class, 4mm, -42m, 422, 4/mmm
            true
        else  # Tetragonal (II) class, 4, -4, 4/m
            c[1, 6] == -c[2, 6]
        end,
    )
end
function symmetry_criteria(::Trigonal, c::StiffnessMatrix)
    return (
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        2c[6, 6] == c[1, 1] - c[1, 2],
        c[1, 4] == -c[2, 4] == -c[5, 6],
        if iszero(c[1, 5])  # Rhombohedral (I) class, 32, -3m, 3m
            true
        else  # Rhombohedral (II) class, 3, -3
            -c[1, 5] == c[2, 5] == c[4, 6]
        end,
    )
end
function symmetry_criteria(::Orthorhombic, c::StiffnessMatrix)
    return (
        all(iszero, (c[1:3, 4:6]..., c[4, 5:6]..., c[5, 6])),
        all(!iszero, (c[1, 1:3]..., c[2, 2:3]..., c[3, 3], c[4, 4], c[5, 5], c[6, 6])),
    )
end
function symmetry_criteria(::Monoclinic, c::StiffnessMatrix)
    return (
        all(iszero, (c[1:3, 4]..., c[5, 6])),
        all(!iszero, (c[1, 1:3]..., c[2, 2:3]..., c[3, 3], c[4, 4], c[5, 5], c[6, 6])),
        if iszero(c[4, 5])  # Diad // x2, standard orientation
            all(iszero, c[1:3, 6]) && all(!iszero, c[1:3, 5]) && !iszero(c[4, 6])
        else
            all(!iszero, c[1:3, 6]) && all(iszero, c[1:3, 5]) && iszero(c[4, 6])
        end,
    )
end
symmetry_criteria(C::CrystalSystem, s::ComplianceMatrix) = symmetry_criteria(C, inv(s))

issystem(x::Union{StiffnessMatrix,ComplianceMatrix}, system::CrystalSystem) =
    all(symmetry_criteria(system, x))

function whichsystem(x::Union{StiffnessMatrix,ComplianceMatrix})
    for system in
        (Cubic(), Hexagonal(), Tetragonal(), Trigonal(), Orthorhombic(), Monoclinic())
        if issystem(x, system)
            return system
        end
    end
    return Triclinic()
end

end
