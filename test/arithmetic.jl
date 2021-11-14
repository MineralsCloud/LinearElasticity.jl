@testset "Test `-`" begin
    @test typeof(-EngineeringStrain(1:6)) == EngineeringStrain{Int}
    @test typeof(-EngineeringStrain(1.0:6.0)) == EngineeringStrain{Float64}
    @test -EngineeringStrain(1:6) == EngineeringStrain(-1, -2, -3, -4, -5, -6)
    @test -EngineeringStress(1.0:6.0) == EngineeringStress(-1.0, -2, -3, -4, -5, -6)
    @test -(-EngineeringStress(1.0:6.0)) == EngineeringStress(1.0:6.0)
end
