begin
    using Revise
    using PolaronQMC
    using Statistics
    using Plots
    using PolaronMobility
    using LaTeXStrings
    using Distributions
end


@time begin
    n_steps = 30000;
    pot = "Frohlich";
    T = 0.2; version=rand(1:10000)
    data_set, energy_arr, error_arr, n_beads = 
        generalPIMC(T, 1.0, 2.0, 3, n_steps, version=1, pot=pot, threads=true)

    println("-----Simulation Ended-----")
    #save("data_arr/MultiThread/$(pot)/T$(T)_nsteps$(n_steps)_v$(version)_beads$(n_beads).jld", "data", data_set)
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