using LinearElasticity
using Test

@testset "LinearElasticity.jl" begin
    include("tensors.jl")
    include("solve.jl")
    include("ULICS.jl")
end
