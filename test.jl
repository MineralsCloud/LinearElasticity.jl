using CrystallographyBase,
    Crystallography,
    LinearElasticity,
    QuantumESPRESSOBase,
    QuantumESPRESSOBase.PWscf,
    QuantumESPRESSOParser.PWscf,
    QuantumESPRESSOFormatter.PWscf
using LinearElasticity.Symmetry: Triclinic, Monoclinic
using LinearElasticity.ULICS
using LinearElasticity.Solve
using Mustache

lattice = Lattice(cell)
strains = ulics(Monoclinic())
lattices = map(strains) do strain
    distortby(lattice, strain)
end
# cells = map(strains) do strain
#     Cell(distortby(lattice, strain), cell.positions, cell.atoms)
# end
map(enumerate(lattices)) do (i, ls)
    open("relax-strain-$i.in", "w") do io
        render(io, str; cell=basisvectors(ls))
    end
end
# map(enumerate(cells)) do (i, cell)
#     open("vcrealx-$i.in", "w") do io
#         render(io, str; ecut=160, cell=basisvectors(Lattice(cell)), atoms=cell.positions, hasecut=true)
#     end
# end

str = read("relax.out", String)
old_lattice = Lattice(parse(CellParametersCard, read("relax.in", String)))
# relaxed_atoms = parseall(AtomicPositionsCard, str)[end]
s0 = parse_stress(str)[2][end] #* u"Ry/bohr^3"
s1_6 = map(1:6) do i
    str = read("realx-strain-$i.out", String)
    parse_stress(str)[2][end] #* u"Ry/bohr^3"
end
Δs1_6 = TensorStress.(s1_6 .- Ref(s0))
lattices16 = map(1:6) do i
    Lattice(parse(CellParametersCard, read("realx-strain-$i.in", String)))
end
qestrains = map(lattices16) do new_lattice
    strainstate(old_lattice, new_lattice)
end
# ap1_6 = map(1:6) do i
#     str = read("realx-strain-$i.out", String)
#     parseall(AtomicPositionsCard, str)[end]
# end
u"GPa".(
    -(StiffnessMatrix(solve_elastic_constants(qestrains, Δs1_6, Triclinic()))) .*
    u"Ry/bohr^3"
)
