begin
    using Revise
    using PolaronQMC
    using Statistics
    using Plots
    using PolaronMobility
    using LaTeXStrings
    include("GeneralPIMC.jl")
end
using Distributions

@time begin
    pot = "Holstein";
    T = 0.05; version=rand(1:10000)
    energy, errors, data_set = 
        general_Holstein_PIMC(T, 1.0, 2.0, 3, 10000000, version=1)

    println("-----Simulation Ended-----")
end