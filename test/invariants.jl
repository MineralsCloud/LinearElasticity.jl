using LinearAlgebra: diagm, norm

# Example from https://www.continuummechanics.org/hydrodeviatoricstrain.html
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

@testset "Test unitary transformation" begin
    E = TensorStrain([
        0.5 0.3 0.2
        0.3 -0.2 -0.1
        0.2 -0.1 0.1
    ])
    evecs = principal_axes(E)
    normal_strains = transpose(evecs) * E * evecs
    @test norm(normal_strains - [
        -0.370577 0 0
        0 0.115308 0
        0 0 0.655269
    ]) < 1e-6
    @test isapprox(principal_values(E), [-0.370577, 0.115308, 0.655269]; atol = 1e-6)
    @test tr(normal_strains) ≈ tr(E) ≈ 2 / 5
end

# Example from https://www.continuummechanics.org/principalstrain.html
@testset "Test principal invariants" begin
    E = TensorStrain([
        0.5 0.3 0.2
        0.3 -0.2 -0.1
        0.2 -0.1 0.1
    ])
    @test length(principal_invariants(E)) == 3
    @test principal_invariants(E)[1] == 0.4
    @test principal_invariants(E)[:I2] ≈ -0.21
    @test principal_invariants(E)[:I3] == -0.028
end
