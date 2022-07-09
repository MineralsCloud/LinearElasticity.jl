using LinearAlgebra: dot, norm
using LinearElasticity.ULICS: ϵ₁, ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆

@testset "Test linear-independency of the strains" begin
    θ(ϵᵢ, ϵⱼ) = rad2deg(acos(dot(ϵᵢ, ϵⱼ) / norm(ϵᵢ) / norm(ϵⱼ)))
    for ϵⱼ in (ϵ₂, ϵ₃, ϵ₄, ϵ₅, ϵ₆)
        @test θ(ϵ₁, ϵⱼ) == 90
    end
    @test θ(ϵ₂, ϵ₃) == 90.62963662344143
    @test θ(ϵ₂, ϵ₄) == 95.67589438805236
    @test θ(ϵ₂, ϵ₅) == 68.73774589826085
    @test θ(ϵ₂, ϵ₆) == 113.99101680223384
    @test θ(ϵ₃, ϵ₄) == 86.85029408329136
    @test θ(ϵ₃, ϵ₅) == θ(ϵ₃, ϵ₆) == 98.21321070173819
    @test θ(ϵ₄, ϵ₅) == 102.0515234011519
    @test θ(ϵ₄, ϵ₆) == 80.51234100434232
    @test θ(ϵ₅, ϵ₆) == 115.37693352515231
end
