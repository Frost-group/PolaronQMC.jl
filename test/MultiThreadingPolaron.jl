begin
    using Revise
    using PolaronQMC
    using Statistics
    using Plots
    using PolaronMobility
    using LaTeXStrings
    include("GeneralPIMC.jl")
end

@time begin
    n_steps = 15000;
    a, b, c = generalPIMC(
                0.2, #Temperature
                1.0, # mass
                1.0, # ω (has to be float)
                2.0, # α (has to be float)
                1, # no of particles
                3, # number of n_dimensions
                "Primitive", # regime type
                true, # fixing beads or not
                0.03, # fixed_τ
                300, # n_beads (if τ is not fixed)
                n_steps, # No. of steps
                50000, # number of thermalisation steps
                "Single", # movers
                "Frohlich", # potential type
                "Virial", # estimators
                false, # Not quick steps
                true, # threading
                1.0, # Start Range
                1,
                1, 
                0.01, # observable skips
                0.2, # equilibrium skips
                1,
                true,
                16
                );
end