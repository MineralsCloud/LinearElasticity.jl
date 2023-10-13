module EffectiveModuli

abstract type Model end
struct Voigt <: Model end
struct Reuss <: Model end
struct VoigtReussHill <: Model end

bulk_modulus(c::StiffnessMatrix, ::Voigt) =
    (c[1, 1] + c[2, 2] + c[3, 3] + 2(c[1, 2] + c[2, 3] + c[3, 1])) / 9
bulk_modulus(s::ComplianceMatrix, ::Voigt) = bulk_modulus(inv(s), Voigt())
bulk_modulus(s::ComplianceMatrix, ::Reuss) =
    inv(s[1, 1] + s[2, 2] + s[3, 3] + 2(s[1, 2] + s[2, 3] + s[3, 1]))
bulk_modulus(c::StiffnessMatrix, ::Reuss) = bulk_modulus(inv(c), Voigt())
bulk_modulus(c::StiffnessMatrix, ::VoigtReussHill) =
    (bulk_modulus(c, Voigt()) + bulk_modulus(inv(c), Reuss())) / 2
bulk_modulus(s::ComplianceMatrix, ::VoigtReussHill) =
    (bulk_modulus(inv(s), Voigt()) + bulk_modulus(s, Reuss())) / 2

shear_modulus(c::StiffnessMatrix, ::Voigt) =
    (
        c[1, 1] + c[2, 2] + c[3, 3] - (c[1, 2] + c[2, 3] + c[3, 1]) +
        3(c[4, 4] + c[5, 5] + c[6, 6])
    ) / 15
shear_modulus(s::ComplianceMatrix, ::Voigt) = shear_modulus(inv(s), Voigt())
shear_modulus(s::ComplianceMatrix, ::Reuss) =
    inv(
        4(s[1, 1] + s[2, 2] + s[3, 3] - (s[1, 2] + s[2, 3] + s[3, 1])) +
        3(s[4, 4] + s[5, 5] + s[6, 6]),
    ) * 15
shear_modulus(c::StiffnessMatrix, ::Reuss) = shear_modulus(inv(c), Voigt())
shear_modulus(c::StiffnessMatrix, ::VoigtReussHill) =
    (shear_modulus(c, Voigt()) + shear_modulus(inv(c), Reuss())) / 2
shear_modulus(s::ComplianceMatrix, ::VoigtReussHill) =
    (shear_modulus(inv(s), Voigt()) + shear_modulus(s, Reuss())) / 2

end
