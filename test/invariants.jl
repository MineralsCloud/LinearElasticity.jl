using LinearAlgebra: diagm

@testset "Test `hydrostatic` and `deviatoric`" begin
    E = TensorStrain([
        0.5 0.3 0.2
        0.3 -0.2 -0.1
        0.2 -0.1 0.1
    ])
    @test hydrostatic(E) == TensorStrain(diagm([2, 2, 2] / 15))
    @test deviatoric(E) == TensorStrain(E - diagm([2, 2, 2] / 15))
    @test TensorStrain(hydrostatic(E) + deviatoric(E)) ≈ E
end
