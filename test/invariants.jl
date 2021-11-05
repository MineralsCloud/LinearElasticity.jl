using LinearAlgebra: I, diagm, norm, tr
using Unitful: @u_str

# Example from https://www.continuummechanics.org/hydrodeviatoricstrain.html
@testset "Test `hydrostatic` and `deviatoric`" begin
    E = TensorStrain([
        0.5 0.3 0.2
        0.3 -0.2 -0.1
        0.2 -0.1 0.1
    ])
    @test hydrostatic(E) == TensorStrain(diagm(0 => [2, 2, 2] / 15))
    @test deviatoric(E) == TensorStrain(E - diagm(0 => [2, 2, 2] / 15))
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
    @test principal_invariants(E)[2] ≈ -0.21
    @test principal_invariants(E)[3] == -0.028
end

@testset "Tests from homework 1" begin
    σ = TensorStress([
        300 100 0
        100 100 0
        0 0 0
    ] * u"MPa")
    @test length(principal_invariants(σ)) == 3
    @test principal_invariants(σ) == (400u"MPa", 2e4u"MPa^2", 0u"MPa^3")
    @test maximum(principal_values(σ)) ≈ (200 + 100 * sqrt(2))u"MPa"
    @test sort(principal_values(σ))[2] ≈ (200 - 100 * sqrt(2))u"MPa"
    @test hydrostatic(σ) == TensorStress(diagm(0 => [400, 400, 400] * u"MPa" / 3))
    @test hydrostatic(σ) + deviatoric(σ) == σ
end

@testset "Tests from homework 2" begin
    σ₂₂ = 200u"MPa"
    σ₁₂ = σ₂₃ = 141u"MPa"
    σ = TensorStress([
        0u"MPa" σ₁₂ 0u"MPa"
        σ₁₂ σ₂₂ σ₂₃
        0u"MPa" σ₂₃ 0u"MPa"
    ])
    @test minimum(principal_values(σ)) ≈ -123.07397876041034u"MPa"
    @test isapprox(sort(principal_values(σ))[2], 0u"MPa"; atol = 1e-12u"MPa")
    @test maximum(principal_values(σ)) == 323.07397876041034u"MPa"
    @test norm(
        principal_axes(σ) - [
            0.601723 0.707107 0.371389
            -0.525223 0 0.850965
            0.601723 -0.707107 0.371389
        ],
    ) < 1e-6
end
