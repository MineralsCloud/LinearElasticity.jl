using LinearAlgebra: dot, norm
using LinearElasticity.ULICS: U₁, U₂, U₃, U₄, U₅, U₆

@testset "Test linear-independency of the strains" begin
    θ(Uᵢ, Uⱼ) = rad2deg(acos(dot(Uᵢ, Uⱼ) / norm(Uᵢ) / norm(Uⱼ)))
    for Uⱼ in (U₂, U₃, U₄, U₅, U₆)
        @test θ(U₁, Uⱼ) == 90
    end
    @test θ(U₂, U₃) == 90.62963662344143
    @test θ(U₂, U₄) == 95.67589438805236
    @test θ(U₂, U₅) == 68.73774589826085
    @test θ(U₂, U₆) == 113.99101680223384
    @test θ(U₃, U₄) == 86.85029408329136
    @test θ(U₃, U₅) == θ(U₃, U₆) == 98.21321070173819
    @test θ(U₄, U₅) == 102.0515234011519
    @test θ(U₄, U₆) == 80.51234100434232
    @test θ(U₅, U₆) == 115.37693352515231
end
