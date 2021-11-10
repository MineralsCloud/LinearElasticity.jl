using CrystallographyBase: Lattice
using LinearAlgebra: I

export distort

# See https://link.springer.com/content/pdf/10.1007%2F978-3-7091-0382-1_7.pdf and https://doi.org/10.2138/am-1997-1-207
distort(lattice::Lattice, strain::TensorStrain) = Lattice((I + strain.data) * lattice.data)
distort(lattice::Lattice, strain::EngineeringStrain) =
    distort(lattice, TensorStrain(strain))
