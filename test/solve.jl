using Crystallography: Hexagonal
using LinearElasticity.SymmetryCriteria: ishexagonal

@testset "Test solving elastic constants on GaN (P6₃mc structure)" begin
    positive_strains = map(1:6) do i
        EngineeringStrain([j == i ? 0.005 : 0 for j in 1:6])
    end
    negative_strains = -positive_strains
    strains = collect(Iterators.flatten(zip(positive_strains, negative_strains)))  # Combine two vectors with alternating strains
    stresses = [
        EngineeringStress(-7.445e-5, -3.1e-7, 1.012e-5, 0, 0, 0),
        EngineeringStress(0.00014751, 7.824e-5, 6.614e-5, 0, 0, 0),
        EngineeringStress(-3.09e-6, -7.17e-5, 1.021e-5, 0, 0, 0),
        EngineeringStress(7.545e-5, 0.00015034, 6.606e-5, 0, 0, 0),
        EngineeringStress(9.66e-6, 9.66e-6, -8.372e-5, 0, 0, 0),
        EngineeringStress(6.551e-5, 6.551e-5, 0.00016129, 0, 0, 0),
        EngineeringStress(3.536e-5, 3.82e-5, 3.79e-5, -3.041e-5, 0, 0),
        EngineeringStress(3.535e-5, 3.821e-5, 3.79e-5, 3.042e-5, 0, 0),
        EngineeringStress(3.512e-5, 3.845e-5, 3.79e-5, 0, -3.039e-5, 0),
        EngineeringStress(3.511e-5, 3.845e-5, 3.79e-5, 0, 3.038e-5, 0),
        EngineeringStress(3.514e-5, 3.862e-5, 3.773e-5, 0, 0, -3.313e-5),
        EngineeringStress(3.51e-5, 3.865e-5, 3.765e-5, 0, 0, 3.858e-5),
    ]
    stiffness_matrix = -solve_elastic_constants(Hexagonal(), strains, stresses)
    @test stiffness_matrix == StiffnessMatrix(
        [
            0.022199649999999998 0.00785485 0.005589249999999999 0 0 0
            0.00785485 0.022199649999999998 0.005589249999999999 0 0 0
            0.005589249999999999 0.005589249999999999 0.02450099999999999 0 0 0
            0 0 0 0.00608 0 0
            0 0 0 0 0.00608 0
            0 0 0 0 0 0.007172399999999999
        ],
    )
    # Reference values: https://materialsproject.org/materials/mp-804?formula=GaN#elastic_constants
    @test stiffness_matrix * u"Ry/bohr^3" .|> u"GPa" |> ustrip == [
        326.56812555486005 115.54883257234204 82.22070599119812 0 0 0
        115.54883257234204 326.56812555486005 82.22070599119812 0 0 0
        82.22070599119812 82.22070599119812 360.42215279158114 0 0 0
        0 0 0 89.43988771775904 0 0
        0 0 0 0 89.43988771775904 0
        0 0 0 0 0 105.509646491259
    ]
    # Reference values: https://materialsproject.org/materials/mp-804?formula=GaN#elastic_constants
    @test inv(stiffness_matrix) * u"bohr^3/Ry" .|> u"TPa^(-1)" |> ustrip == [
        3.6052277063145053 -1.1336754805274927 -0.5638187534378214 0 0 0
        -1.1336754805274927 3.6052277063145053 -0.5638187534378214 0 0 0
        -0.5638187534378214 -0.5638187534378214 3.031764677764818 0 0 0
        0 0 0 11.180693821482086 0 0
        0 0 0 0 11.180693821482086 0
        0 0 0 0 0 9.477806373683997
    ]
    @test ishexagonal(stiffness_matrix)
end

@testset "Test solving elastic constants on KNO₂ (Cm structure)" begin
    positive_strains = map(1:6) do i
        EngineeringStrain([j == i ? 0.005 : 0 for j in 1:6])
    end
    negative_strains = -positive_strains
    strains = collect(Iterators.flatten(zip(positive_strains, negative_strains)))
    stresses = [
        EngineeringStress(-2.07e-6, 1.019e-5, 6.84e-6, 0, -6.3e-7, 0),
        EngineeringStress(2.545e-5, 1.552e-5, 1.078e-5, 0, 2.4e-7, 0),
        EngineeringStress(8.61e-6, 4.6e-7, 3.38e-6, 0, -1.1e-7, 0),
        EngineeringStress(1.386e-5, 2.551e-5, 1.456e-5, 0, -3.0e-7, 0),
        EngineeringStress(9.27e-6, 7.37e-6, 1.6e-6, 0, -9.0e-8, 0),
        EngineeringStress(1.322e-5, 1.857e-5, 1.621e-5, 0, -3.5e-7, 0),
        EngineeringStress(1.12e-5, 1.281e-5, 8.8e-6, -2.17e-6, -1.9e-7, -3.0e-8),
        EngineeringStress(1.121e-5, 1.282e-5, 8.8e-6, 2.17e-6, -1.9e-7, 3.0e-8),
        EngineeringStress(1.076e-5, 1.3e-5, 8.97e-6, 0, -6.1e-7, 0),
        EngineeringStress(1.171e-5, 1.271e-5, 8.61e-6, 0, 2.3e-7, 0),
        EngineeringStress(1.121e-5, 1.292e-5, 8.85e-6, 5.0e-8, -2.1e-7, -2.26e-6),
        EngineeringStress(1.12e-5, 1.292e-5, 8.85e-6, -5.0e-8, -2.1e-7, 2.26e-6),
    ]
    # Reference values: https://materialsproject.org/materials/mp-34857?formula=KNO2#elastic_constants
end
