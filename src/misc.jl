using Crystallography: Lattice
using LinearAlgebra: I

export distortby, strainstate, isuniaxial, isbiaxial

# See https://link.springer.com/content/pdf/10.1007%2F978-3-7091-0382-1_7.pdf and https://doi.org/10.2138/am-1997-1-207
distortby(lattice::Lattice, strain::TensorStrain) =
    Lattice((I + strain.data) * lattice.data)
distortby(lattice::Lattice, strain::EngineeringStrain) =
    distortby(lattice, TensorStrain(strain))

strainstate(old::Lattice, new::Lattice) = TensorStrain(new.data / old.data - I)

# No shear components, only one normal component is nonzero
isuniaxial(x::Union{EngineeringStress,EngineeringStrain}) =
    iszero(x[4:end]) && length(filter(iszero, x[1:3])) == 2
isuniaxial(σ::TensorStress) = isuniaxial(EngineeringStress(σ))
isuniaxial(ε::TensorStrain) = isuniaxial(EngineeringStrain(ε))

function isbiaxial(x::Union{EngineeringStress,EngineeringStrain})
    n = length(filter(!iszero, x[4:end]))
    return if n > 1
        false  # Triaxial
    elseif n == 1
        if all(iszero, x[1:3])  # Pure shear
            true
        else
            if !iszero(x[6])  # 12 ≠ 0
                !iszero(x[1]) && !iszero(x[2]) && iszero(x[3])
            elseif !iszero(x[5])  # 13 ≠ 0
                !iszero(x[1]) && iszero(x[2]) && !iszero(x[3])
            else  # 23 ≠ 0
                iszero(x[1]) && !iszero(x[2]) && !iszero(x[3])
            end
        end
    else  # `n` is zero, no shear components
        length(filter(iszero, x[1:3])) == 1  # Only two normal components
    end
end
isbiaxial(σ::TensorStress) = isbiaxial(EngineeringStress(σ))
isbiaxial(ε::TensorStrain) = isbiaxial(EngineeringStrain(ε))
