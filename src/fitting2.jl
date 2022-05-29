using Crystallography: CrystalSystem, Cubic, Hexagonal, Tetragonal, Orthorhombic
using LinearAlgebra: Symmetric, dot

export fit_elastic_constant

function fit_elastic_constant(
    ::Cubic,
    strains::AbstractVector{<:EngineeringStrain},
    stresses::AbstractVector{<:EngineeringStress},
)

end
