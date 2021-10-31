using Crystallography:
    CrystalSystem, Cubic, Hexagonal, Tetragonal, Trigonal, Orthorhombic, Monoclinic
using LinearAlgebra: issymmetric

export issystem

function symmetry_criteria(::Cubic, c::StiffnessMatrix)
    return all((
        issymmetric(c),
        c[1, 1] == c[2, 2] == c[3, 3],
        c[4, 4] == c[5, 5] == c[6, 6],
        c[1, 2] == c[1, 3] == c[2, 3],
    ))
end
function symmetry_criteria(::Hexagonal, c::StiffnessMatrix)
    return all((
        issymmetric(c),
        c[1, 1] == c[2, 2],
        c[4, 4] == c[5, 5],
        c[1, 3] == c[2, 3],
        2c[6, 6] == c[1, 1] - c[1, 2],
    ))
end
function symmetry_criteria(::Tetragonal, c::StiffnessMatrix)
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
function symmetry_criteria(::Trigonal, c::StiffnessMatrix)
    return all((
        issymmetric(c),
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
    ))
end
function symmetry_criteria(::Orthorhombic, c::StiffnessMatrix)
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
function symmetry_criteria(::Monoclinic, c::StiffnessMatrix)
    return all((
        issymmetric(c),
        all(
            iszero,
            (c[1, 4], c[2, 4], c[3, 4], c[4, 5], c[1, 6], c[2, 6], c[3, 6], c[5, 6]),
        ),
    ))
end
symmetry_criteria(C::CrystalSystem, s::ComplianceMatrix) = symmetry_criteria(C, inv(s))

issystem(C::CrystalSystem, x::Union{StiffnessMatrix,ComplianceMatrix}) =
    all(symmetry_criteria(C, x))
