@testset "Test if `isuniaxial` works" begin
    @testset "Simple tests" begin
        @test isuniaxial(EngineeringStrain(0, 1e-5, 0, 0, 0, 0))
        @test !isuniaxial(EngineeringStrain(0, 0, 0, 0, 1e-5, 0))
        @test !isuniaxial(EngineeringStrain(1e-5, 1e-5, 0, 0, 0, 0))
        @test !isuniaxial(EngineeringStress(0, 1e-5, 2e-5, 0, 0, 0))
        @test !isuniaxial(EngineeringStress(0, 1e-5, 0, 0, 2e-6, 0))
    end
    σ₂₂ = 200
    σ₁₂ = σ₂₃ = 141
    @test isuniaxial(TensorStress([
        0 0 0
        0 σ₂₂ 0
        0 0 0
    ]))
    @test !isuniaxial(TensorStress([
        0 σ₁₂ 0
        σ₁₂ σ₂₂ 0
        0 0 0
    ]))
    @test !isuniaxial(TensorStress([
        0 0 0
        0 0 σ₂₃
        0 σ₂₃ 0
    ]))
    σ = TensorStress([
        0 σ₁₂ 0
        σ₁₂ σ₂₂ σ₂₃
        0 σ₂₃ 0
    ])
    @test !isuniaxial(σ)
end

# See https://www.idc-online.com/technical_references/pdfs/mechanical_engineering/Plane_Stress_Tensor_Biaxial_Stress.pdf
@testset "Test if `isbiaxial` works" begin
    @testset "Simple tests" begin
        @test !isbiaxial(EngineeringStrain(0, 1e-5, 0, 0, 0, 0))
        @test isbiaxial(EngineeringStrain(0, 0, 0, 0, 1e-5, 0))
        @test isbiaxial(EngineeringStrain(1e-5, 1e-5, 0, 0, 0, 0))
        @test isbiaxial(EngineeringStress(0, 1e-5, 2e-5, 0, 0, 0))
        @test !isbiaxial(EngineeringStress(0, 1e-5, 0, 0, 2e-6, 0))
        @test !isbiaxial(EngineeringStress(1e-5, 0, 1e-5, 0, 2e-6, 0))
    end
    σ₂₂ = 200
    σ₁₂ = σ₂₃ = 141
    @test isbiaxial(TensorStress([
        σ₂₂/2 0 0
        0 σ₂₂ 0
        0 0 0
    ]))  # Without shear
    @test isbiaxial(TensorStress([
        σ₂₂/2 σ₁₂ 0
        σ₁₂ σ₂₂ 0
        0 0 0
    ]))  # With shear
    @test !isbiaxial(TensorStress([
        σ₂₂/2 σ₁₂ 0
        σ₁₂ σ₂₂ σ₂₃
        0 σ₂₃ 0
    ]))  # With shear
    @test !isbiaxial(TensorStress([
        σ₂₂/2 0 0
        0 σ₂₂ σ₂₃
        0 σ₂₃ 0
    ]))  # With shear
    @test isbiaxial(TensorStress([
        0 0 0
        0 σ₂₂ σ₂₃
        0 σ₂₃ σ₂₂/2
    ]))  # With shear
    @test isbiaxial(TensorStress([
        0 0 0
        0 0 σ₂₃
        0 σ₂₃ 0
    ]))  # Pure shear
    @test !isbiaxial(TensorStress([
        0 σ₁₂ 0
        σ₁₂ 0 σ₂₃
        0 σ₂₃ 0
    ]))  # Pure shear
    @test isbiaxial(TensorStress([
        0 σ₁₂ 0
        σ₁₂ 0 0
        0 0 0
    ]))  # Pure shear
    σ = TensorStress([
        0 σ₁₂ 0
        σ₁₂ σ₂₂ σ₂₃
        0 σ₂₃ 0
    ])
    @test !isbiaxial(σ)
    @test !isbiaxial(EngineeringStress(0, σ₂₂, 0, σ₂₃, 0, σ₁₂))
end
