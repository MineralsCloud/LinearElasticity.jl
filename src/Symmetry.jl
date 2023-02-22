module Symmetry

using LinearElasticityBase: StiffnessMatrix, ComplianceMatrix

export Cubic, Hexagonal, Tetragonal, Trigonal, Orthorhombic, Monoclinic, Triclinic
export hassymmetry, guesssymmetry, isisotropic

abstract type SymmetryConstraint end
struct Triclinic <: SymmetryConstraint end
struct Monoclinic <: SymmetryConstraint end
struct Orthorhombic <: SymmetryConstraint end
struct Tetragonal <: SymmetryConstraint end
struct Cubic <: SymmetryConstraint end
struct Trigonal <: SymmetryConstraint end
struct Hexagonal <: SymmetryConstraint end

function meetcriteria(x::Union{StiffnessMatrix,ComplianceMatrix}, ::Cubic)
    return (
        all(iszero, x[1:3, 4:6]),
        all(iszero, x[4, 5:6]),
        iszero(x[5, 6]),
        all(!iszero, (x[1, 1], x[4, 4], x[1, 3])),
        x[1, 1] == x[2, 2] == x[3, 3],
        x[4, 4] == x[5, 5] == x[6, 6],
        x[1, 2] == x[1, 3] == x[2, 3],
    )
end
function meetcriteria(c::StiffnessMatrix, ::Hexagonal)
    return (
        all(iszero, c[1:3, 4:6]),
        all(iszero, c[4, 5:6]),
        iszero(c[5, 6]),
        all(!iszero, (c[1, 1], c[3, 3], c[4, 4], c[1, 2], c[1, 3], c[6, 6])),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        2c[6, 6] == c[1, 1] - c[1, 2],
    )
end
function meetcriteria(s::ComplianceMatrix, ::Hexagonal)
    return (
        all(iszero, s[1:3, 4:6]),
        all(iszero, s[4, 5:6]),
        iszero(s[5, 6]),
        all(!iszero, (s[1, 1], s[3, 3], s[4, 4], s[1, 2], s[1, 3], s[6, 6])),
        s[1, 1] == s[2, 2],
        s[4, 4] == s[5, 5],
        s[1, 3] == s[2, 3],
        s[6, 6] == 2(s[1, 1] - s[1, 2]),
    )
end
function meetcriteria(x::Union{StiffnessMatrix,ComplianceMatrix}, ::Tetragonal)
    return (
        all(iszero, x[1:3, 4:5]),
        all(iszero, x[4, 5:6]),
        all(iszero, (x[3, 6], x[5, 6])),
        all(!iszero, (x[1, 1], x[1, 2], x[1, 3], x[3, 3], x[4, 4], x[6, 6])),
        x[1, 1] == x[2, 2],
        x[4, 4] == x[5, 5],
        x[1, 3] == x[2, 3],
        if all(iszero, (x[2, 6], x[1, 6]))  # Tetragonal (I) class, 4mm, -42m, 422, 4/mmm
            true
        else  # Tetragonal (II) class, 4, -4, 4/m
            x[1, 6] == -x[2, 6]
        end,
    )
end
function meetcriteria(c::StiffnessMatrix, ::Trigonal)
    return (
        all(iszero, c[1:3, 6]),
        all(iszero, c[3, 4:5]),
        iszero(c[4, 5]),
        all(!iszero, c[1, 1:4]),
        all(!iszero, c[2, 2:3]),
        all(!iszero, (c[3, 3], c[4, 4], c[6, 6])),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        2c[6, 6] == c[1, 1] - c[1, 2],
        c[1, 4] == -c[2, 4] == c[5, 6],
        if iszero(c[1, 5])  # Rhombohedral (I) class, 32, -3m, 3m
            all(iszero, (c[2, 5], c[4, 6]))
        else  # Rhombohedral (II) class, 3, -3
            -c[1, 5] == c[2, 5] == c[4, 6]
        end,
    )
end
function meetcriteria(s::ComplianceMatrix, ::Trigonal)
    return (
        all(iszero, s[1:3, 6]),
        all(iszero, s[3, 4:5]),
        iszero(s[4, 5]),
        all(!iszero, s[1, 1:4]),
        all(!iszero, s[2, 2:3]),
        all(!iszero, (s[3, 3], s[4, 4], s[6, 6])),
        s[1, 1] == s[2, 2],
        s[4, 4] == s[5, 5],
        s[1, 3] == s[2, 3],
        s[6, 6] == 2(s[1, 1] - s[1, 2]),
        2s[1, 4] == -2s[2, 4] == s[5, 6],
        if iszero(s[1, 5])  # Rhombohedral (I) class, 32, -3m, 3m
            all(iszero, (s[2, 5], s[4, 6]))
        else  # Rhombohedral (II) class, 3, -3
            -2s[1, 5] == 2s[2, 5] == s[4, 6]
        end,
    )
end
function meetcriteria(x::Union{StiffnessMatrix,ComplianceMatrix}, ::Orthorhombic)
    return (
        all(iszero, x[1:3, 4:6]),
        all(iszero, x[4, 5:6]),
        iszero(x[5, 6]),
        all(!iszero, x[1:2, 1:3]),
        all(!iszero, (x[3, 3], x[4, 4], x[5, 5], x[6, 6])),
    )
end
function meetcriteria(x::Union{StiffnessMatrix,ComplianceMatrix}, ::Monoclinic)
    return (
        all(iszero, x[1:3, 4]),
        iszero(x[5, 6]),
        all(!iszero, x[1:2, 1:3]),
        all(!iszero, (x[3, 3], x[4, 4], x[5, 5], x[6, 6])),
        if iszero(x[4, 5])  # Diad // x2, standard orientation
            all(iszero, x[1:3, 6]) && all(!iszero, x[1:3, 5]) && !iszero(x[4, 6])
        else
            all(!iszero, x[1:3, 6]) && all(iszero, x[1:3, 5]) && iszero(x[4, 6])
        end,
    )
end
meetcriteria(x::Union{StiffnessMatrix,ComplianceMatrix}, ::Triclinic) =
    all(!iszero, x.data.data)

hassymmetry(x::Union{StiffnessMatrix,ComplianceMatrix}, cstr::SymmetryConstraint) =
    all(meetcriteria(cstr, x))

function guesssymmetry(x::Union{StiffnessMatrix,ComplianceMatrix})
    for symmetry in (
        Cubic(),
        Hexagonal(),
        Tetragonal(),
        Trigonal(),
        Orthorhombic(),
        Monoclinic(),
        Triclinic(),
    )
        if hassymmetry(x, symmetry)
            return symmetry
        end
    end
end

function isisotropic(c::StiffnessMatrix)
    return all((
        all(iszero, c[1:3, 4:6]),
        all(iszero, c[4, 5:6]),
        iszero(c[5, 6]),
        all(!iszero, (c[1, 1], c[4, 4], c[1, 3])),
        c[1, 1] == c[2, 2] == c[3, 3],
        2c[4, 4] == 2c[5, 5] == 2c[6, 6] == c[1, 1] - c[1, 2],
        c[1, 2] == c[1, 3] == c[2, 3],
    ))
end
function isisotropic(s::ComplianceMatrix)
    return all((
        all(iszero, s[1:3, 4:6]),
        all(iszero, s[4, 5:6]),
        iszero(s[5, 6]),
        all(!iszero, (s[1, 1], s[4, 4], s[1, 3])),
        s[1, 1] == s[2, 2] == s[3, 3],
        s[4, 4] == s[5, 5] == s[6, 6] == 2(s[1, 1] - s[1, 2]),
        s[1, 2] == s[1, 3] == s[2, 3],
    ))
end

end
