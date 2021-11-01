using LinearElasticity
using Test

@testset "LinearElasticity.jl" begin
    include("tensors.jl")
    include("invariants.jl")
end
