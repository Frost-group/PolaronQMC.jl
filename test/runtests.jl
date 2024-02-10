using PolaronQMC

using Test

@testset verbose = true "PolaronQMC" begin
    @testset "Harmonic" begin
        include("harmonic.jl")
    end

    @testset "Frohlich" begin
        include("FrohlichTest.jl")
    end

    @testset "Holstein" begin
        include("HolsteinTest.jl")
    end
end

