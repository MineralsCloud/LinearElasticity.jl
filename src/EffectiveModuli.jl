module EffectiveModuli

using ..LinearElasticity: StiffnessMatrix, ComplianceMatrix

import ..Isotropic: bulk_modulus, shear_modulus

export Voigt, Reuss, VoigtReussHill, bulk_modulus, shear_modulus

abstract type Model end
struct Voigt <: Model end
struct Reuss <: Model end
struct VoigtReussHill <: Model end

"""
    bulk_modulus(c::StiffnessMatrix, ::Voigt)
    bulk_modulus(s::ComplianceMatrix, ::Voigt)

Calculate the bulk modulus of a material from its stiffness matrix by:

```math
K_V = \\frac{c_{11} + c_{22} + c_{33} + 2(c_{12} + c_{23} + c_{31})}{9}.
```
"""
bulk_modulus(c::StiffnessMatrix, ::Voigt) =
    (c[1, 1] + c[2, 2] + c[3, 3] + 2(c[1, 2] + c[2, 3] + c[3, 1])) / 9
bulk_modulus(s::ComplianceMatrix, ::Voigt) = bulk_modulus(inv(s), Voigt())
"""
    bulk_modulus(s::ComplianceMatrix, ::Reuss)
    bulk_modulus(c::StiffnessMatrix, ::Reuss)

Calculate the bulk modulus of a material from its compliance matrix by:

```math
K_R = \\frac{1}{s_{11} + s_{22} + s_{33} + 2(s_{12} + s_{23} + s_{31})}.
```
"""
bulk_modulus(s::ComplianceMatrix, ::Reuss) =
    inv(s[1, 1] + s[2, 2] + s[3, 3] + 2(s[1, 2] + s[2, 3] + s[3, 1]))
bulk_modulus(c::StiffnessMatrix, ::Reuss) = bulk_modulus(inv(c), Voigt())
"""
    bulk_modulus(c::StiffnessMatrix, ::VoigtReussHill)
    bulk_modulus(s::ComplianceMatrix, ::VoigtReussHill)

Calculate the bulk modulus of a material from its stiffness or compliance matrix by:

```math
K_{VRH} = \\frac{K_V + K_R}{2}.
```
"""
bulk_modulus(c::StiffnessMatrix, ::VoigtReussHill) =
    (bulk_modulus(c, Voigt()) + bulk_modulus(inv(c), Reuss())) / 2
bulk_modulus(s::ComplianceMatrix, ::VoigtReussHill) =
    (bulk_modulus(inv(s), Voigt()) + bulk_modulus(s, Reuss())) / 2

"""
    shear_modulus(c::StiffnessMatrix, ::Voigt)

Calculate the shear modulus of a material from its stiffness matrix by:

```math
G_V = \\frac{(c_{11} + c_{22} + c_{33}) - (c_{12} + c_{23} + c_{31}) + 3(c_{44} + c_{55} + c_{66})}{15}.
```
"""
shear_modulus(c::StiffnessMatrix, ::Voigt) =
    (
        c[1, 1] + c[2, 2] + c[3, 3] - (c[1, 2] + c[2, 3] + c[3, 1]) +
        3(c[4, 4] + c[5, 5] + c[6, 6])
    ) / 15
shear_modulus(s::ComplianceMatrix, ::Voigt) = shear_modulus(inv(s), Voigt())
"""
    shear_modulus(s::ComplianceMatrix, ::Reuss)

Calculate the shear modulus of a material from its compliance matrix by:

```math
G_R = \\frac{15}{4(s_{11} + s_{22} + s_{33}) - 4(s_{12} + s_{23} + s_{31}) + 3(s_{44} + s_{55} + s_{66})}.
```
"""
shear_modulus(s::ComplianceMatrix, ::Reuss) =
    inv(
        4(s[1, 1] + s[2, 2] + s[3, 3] - (s[1, 2] + s[2, 3] + s[3, 1])) +
        3(s[4, 4] + s[5, 5] + s[6, 6]),
    ) * 15
shear_modulus(c::StiffnessMatrix, ::Reuss) = shear_modulus(inv(c), Voigt())
"""
    shear_modulus(c::StiffnessMatrix, ::VoigtReussHill)
    shear_modulus(s::ComplianceMatrix, ::VoigtReussHill)

Calculate the shear modulus of a material from its stiffness or compliance matrix by:

```math
G_{VRH} = \\frac{G_V + G_R}{2}.
```
"""
shear_modulus(c::StiffnessMatrix, ::VoigtReussHill) =
    (shear_modulus(c, Voigt()) + shear_modulus(inv(c), Reuss())) / 2
shear_modulus(s::ComplianceMatrix, ::VoigtReussHill) =
    (shear_modulus(inv(s), Voigt()) + shear_modulus(s, Reuss())) / 2

end
