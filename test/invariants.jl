using LinearAlgebra: I, diagm, norm, tr
# using Unitful: @u_str

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
    @test transpose(evecs) * evecs ≈ evecs * transpose(evecs) ≈ I
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

# @testset "Tests from homework 1" begin
#     σ = TensorStress([
#         300 100 0
#         100 100 0
#         0 0 0
#     ] * u"MPa")
#     @test principal_values(E)[1] == (200 + 100 √ 2) * u"MPa"
#     @test principal_values(E)[2] == (200 - 100 √ 2) * u"MPa"
# end

@testset "Tests from homework 2" begin
    σ₂₂ = 200
    σ₁₂ = σ₂₃ = 141
    σ = TensorStress([
        0 σ₁₂ 0
        σ₁₂ σ₂₂ σ₂₃
        0 σ₂₃ 0
    ])
    @test principal_values(σ)[1] == -123.07397876041034
    @test isapprox(principal_values(σ)[2], 0; atol = 1e-13)
    @test principal_values(σ)[3] == 323.07397876041034
    @test norm(principal_axes(σ) - [
        0.601723   0.707107     0.371389
        -0.525223   0  0.850965
         0.601723  -0.707107     0.371389
    ]) < 1e-6
end
