@testset "Test `-`" begin
    @test typeof(-EngineeringStrain(1:6)) == EngineeringStrain{Int}
    @test typeof(-EngineeringStrain(1.0:6.0)) == EngineeringStrain{Float64}
    @test -EngineeringStrain(1:6) == EngineeringStrain(-1, -2, -3, -4, -5, -6)
    @test -EngineeringStress(1.0:6.0) == EngineeringStress(-1.0, -2, -3, -4, -5, -6)
    @test -(-EngineeringStress(1.0:6.0)) == EngineeringStress(1.0:6.0)
end

@testset "Test `*` and `/`" begin
    a = EngineeringStrain(1:6)
    @test typeof(2a) == EngineeringStrain{Int}
    @test typeof(2.0a) == EngineeringStrain{Float64}
    @test 2a == 2 * a == a * 2 == 2.0a == EngineeringStress(2:2:12)
    @test a / 2 == EngineeringStress(0.5, 1, 1.5, 2, 2.5, 3)
    @test 2a / 2 == a
end
