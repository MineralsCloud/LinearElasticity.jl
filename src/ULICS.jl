module ULICS

using ..LinearElasticity: CrystalSystem, EngineeringStrain, minimal_npairs
using Random: shuffle

export ulics

# See Yu, R., Zhu, J. & Ye, H. Q. Calculations of single-crystal elastic constants made simple. Comput Phys Commun 181, 671–675 (2010).
const U₁ = EngineeringStrain(1:6) * 1e-3
const U₂ = EngineeringStrain(-2, 1, 4, -3, 6, -5) * 1e-3
const U₃ = EngineeringStrain(3, -5, -1, 6, 2, -4) * 1e-3
const U₄ = EngineeringStrain(-4, -6, 5, 1, -3, 2) * 1e-3
const U₅ = EngineeringStrain(5, 4, 6, -2, -1, -3) * 1e-3
const U₆ = EngineeringStrain(-6, 3, -2, 5, -4, 1) * 1e-3
const U = [U₁, U₂, U₃, U₄, U₅, U₆]

ulics(system::CrystalSystem) = shuffle(U)[1:minimal_npairs(system)]

end
