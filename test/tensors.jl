using LinearElasticityBase: SymmetricFourthOrderTensor

using LinearElasticity: Cubic
using LinearElasticity.SymmetryCriteria: whichsystem, isisotropic

# Compared with https://ferrite-fem.github.io/Tensors.jl/stable/demos/
@testset "Creating the linear elasticity tensor" begin
    E = 2
    ν = 0.3
    dim = 3
    λ = E * ν / ((1 + ν) * (1 - 2ν))
    μ = E / (2(1 + ν))
    δ(i, j) = i == j ? 1.0 : 0.0
    f = (i, j, k, l) -> λ * δ(i, j) * δ(k, l) + μ * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k))
    C = StiffnessTensor(SymmetricFourthOrderTensor{3}(f))
    @test isapprox(
        C[:, :, 1, 1],
        [
            2.69231 0.0 0.0
            0.0 1.15385 0.0
            0.0 0.0 1.15385
        ];
        rtol=1e-5,
    )
    @test isapprox(
        C[:, :, 2, 3],
        [
            0.0 0.0 0.0
            0.0 0.0 0.769231
            0.0 0.769231 0.0
        ];
        rtol=1e-6,
    )
    c = StiffnessMatrix(C)
    @test isapprox(
        c,
        [
            2.69231 1.15385 1.15385 0.0 0.0 0.0
            1.15385 2.69231 1.15385 0.0 0.0 0.0
            1.15385 1.15385 2.69231 0.0 0.0 0.0
            0.0 0.0 0.0 0.769231 0.0 0.0
            0.0 0.0 0.0 0.0 0.769231 0.0
            0.0 0.0 0.0 0.0 0.0 0.769231
        ];
        rtol=1e-5,
    )  # This is an isotropic system
    @test StiffnessTensor(c) == C
    S = inv(C)
    s = ComplianceMatrix(S)
    @test s ≈ inv(c)
    @test inv(s) ≈ StiffnessMatrix(C)
    @test whichsystem(c) == Cubic()
    @testset "Test for an isotropic system" begin
        @test c[1, 1] * s[1, 1] + 2c[1, 2] * s[1, 2] == 1
        @test c[1, 1] * s[1, 2] + c[1, 2] * s[1, 1] + c[1, 2] * s[1, 2] < eps()
        @test c[4, 4] * s[4, 4] ≈ 1
        @test isisotropic(c)
    end
end

@testset "Creating the stresses and strains" begin
    σ₂₂ = 200
    σ₁₂ = σ₂₃ = 141
    σ = TensorStress([
        0 σ₁₂ 0
        σ₁₂ σ₂₂ σ₂₃
        0 σ₂₃ 0
    ])
    @test EngineeringStress(σ) == EngineeringStress([0, 200, 0, 141, 0, 141])
    @test TensorStress(EngineeringStress(σ)) == σ
    E = TensorStrain([
        0.5 0.3 0.2
        0.3 -0.2 -0.1
        0.2 -0.1 0.1
    ])
    @test EngineeringStrain(E) == EngineeringStrain([0.5, -0.2, 0.1, -0.2, 0.4, 0.6])
    @test TensorStrain(EngineeringStrain(E)) == E
end
