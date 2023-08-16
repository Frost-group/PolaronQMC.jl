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
    n_steps = 100000;
    pot = "Frohlich";
    T = 0.4; version=rand(1:10000)
    data_set, energy_arr, error_arr, n_beads = generalPIMC(
                T, #Temperature
                1.0, # mass
                1.0, # ω (has to be float)
                8.0, # α (has to be float)
                1, # no of particles
                3, # number of n_dimensions
                "Primitive", # regime type
                false, # fixing beads or not
                0.015, # fixed_τ
                200, # n_beads (if τ is not fixed)
                n_steps, # No. of stepsx
                20000, # number of thermalisation steps
                "Single", # movers
                pot, # potential type
                "Virial", # estimators
                false, # Not quick steps
                true, # threading
                1.0, # Start Range
                1,
                1, 
                0.01, # observable skips
                0.5, # equilibrium skips
                1,
                true,
                15
                );

    println("-----Simulation Ended-----")
    save("data_arr/MultiThread/$(pot)/T$(T)_nsteps$(n_steps)_v$(version)_beads$(n_beads).jld", "data", data_set)
end

begin
    particleIndex, dimensionIndex = 1, 1;
    position_all = []
    for i in 1:length(data_set)
        pos_individual = data_set[i]["Position:p$(particleIndex)d$(dimensionIndex)"]
        collect(Iterators.flatten(pos_individual))
        push!(position_all, pos_individual)
    end
    positions_flattened = collect(Iterators.flatten(position_all))
    positions_flattened = collect(Iterators.flatten(positions_flattened))
    positionplot = histogram(positions_flattened, xlab = "Position of runs")
    A2 = Distributions.fit(Normal, positions_flattened)
    println(A2)
    display(positionplot)
end