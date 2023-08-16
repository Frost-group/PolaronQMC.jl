using Test

#using PolaronQMC

function foo()
    return 3
end

println(foo())
#@assert 1==1

#PolaronQMC.Path(nbeads, nparticles,  spatialdimensions=spatialdimensions)


#@testset "BenchmarkMoves.jl" begin include("BenchmarkMoves.jl") end


#@testset "BenchmarkMoves.jl" begin include("BenchmarkMoves.jl") end
a = 3
@test π≈3.14 atol=0.01 #quick one-line test
@test 2 + 2 ≈ 6 atol=1 broken=true
@testset foo()

@testset "trigonometric identities" begin
    θ = 2/3*π
    @test sin(-θ) ≈ -sin(θ)
    @test cos(-θ) ≈ cos(θ)
    @test sin(2θ) ≈ 2*sin(θ)*cos(θ)
    @test cos(2θ) ≈ cos(θ)^2 - sin(θ)^2
end;

@tes