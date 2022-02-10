using Crystallography: Lattice
using LinearAlgebra: I

export distortby, distort, strainstate

# See https://link.springer.com/content/pdf/10.1007%2F978-3-7091-0382-1_7.pdf and https://doi.org/10.2138/am-1997-1-207
distortby(lattice::Lattice, strain::TensorStrain) =
    Lattice((I + strain.data) * lattice.data)
distortby(lattice::Lattice, strain::EngineeringStrain) =
    distortby(lattice, TensorStrain(strain))
const distort = distortby  # For the sake of compatibility

strainstate(old::Lattice, new::Lattice) = TensorStrain(new.data / old.data - I)
