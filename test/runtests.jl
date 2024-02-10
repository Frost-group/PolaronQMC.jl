using PolaronQMC

using Test

@testset verbose = true "PolaronQMC" begin
    @testset "Frohlich" begin
        include("FrohlichTest.jl")
    end
    @testset "Holstein" begin
        include("HolsteinTest.jl")
    end
end

