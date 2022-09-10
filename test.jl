using QuantumESPRESSO.Inputs.PWscf
using QuantumESPRESSO.Outputs.PWscf
using LinearElasticity
using Crystallography

old_cell = read("elast.in", PWInput).cell_parameters |> Lattice

stresses, strains = EngineeringStress[], EngineeringStrain[]
for f in [-0.01, 0.01]
    for x in 1:6
        path = joinpath("e" * string(x), string(f))
        cd(path) do
            stress =
                parse_stress(read("elast.out", String))[2][end] |>
                TensorStress |>
                EngineeringStress
            push!(stresses, stress)
            new_cell = read("elast.in", PWInput).cell_parameters |> Lattice
            strain = strainstate(old_cell, new_cell) |> TensorStrain |> EngineeringStrain
            println(strain)
            push!(strains, strain)
        end
    end
end
-solve_elastic_constants(Tetragonal(), strains, stresses) * u"Ry/bohr^3" .|> u"GPa"
