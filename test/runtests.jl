using PolaronQMC

using Test

@testset "PolaronQMC" begin
    @testset "Frohlich-singlethread" begin
        include("FrohlichSingleThread.jl")
    end
end

