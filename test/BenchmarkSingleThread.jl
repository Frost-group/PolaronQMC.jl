using Revise
using BenchmarkTools
#using PolaronQMC

include("GeneralPIMC.jl")

#=
@time generalPIMC(
    0.9, #Temperature
    1.0, # mass
    1.0, # ω (has to be float)
    9.0, # α (has to be float)
    1, # no of particles
    3, # number of n_dimensions
    "Primitive", # regime type
    true, # fixing beads or not
    0.03, # fixed_τ
    50, # n_beads (if τ is not fixed)
    5000, # No. of steps
    10000, # number of thermalisation steps
    "Single", # movers
    "Frohlich", # potential type
    "Thermodynamic", # estimators
    false, # Not quick steps
    false, # Not threading
    1.0, # Start Range
    1,
    1, 
    0.01, # observable skips
    0.2, # equilibrium skips
    1,
    true,
    1
    )
println("---end---")
=#


@time generalPIMC(
    0.5, #Temperature
    1.0, # mass
    1.0, # ω (has to be float)
    4.0, # α (has to be float)
    1, # no of particles
    3, # number of n_dimensions
    "Primitive", # regime type
    true, # fixing beads or not
    0.03, # fixed_τ
    100, # n_beads (if τ is not fixed)
    1000, # No. of steps
    10000, # number of thermalisation steps
    "Single", # movers
    "Frohlich", # potential type
    "Virial", # estimators
    false, # Not quick steps
    false, # Not threading
    1.0, # Start Range
    1,
    1, 
    0.01, # observable skips
    0.2, # equilibrium skips
    1,
    true,
    1
    )
println("---end---")



