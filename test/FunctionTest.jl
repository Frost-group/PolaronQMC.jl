begin
    using Revise
    using PolaronQMC
    using Statistics
    using Plots
    using PolaronMobility
    using LaTeXStrings
    using BenchmarkTools
    include("GeneralPIMC.jl")
end

@time begin
    n_steps = 100000;
    T = 0.01;
    mean_energy, jacknife_errors, comparison_energy, energies, positions, correlations, acceptance_rates, adjuster_values, equilibrium_skip, observables_skip, version, potential, n_beads, data = generalPIMC(
                T, #Temperature
                1.0, # mass
                1.0, # ω (has to be float)
                7.0, # α (has to be float)
                1, # no of particles
                3, # number of n_dimensions
                "Primitive", # regime type
                true, # fixing beads or not
                0.04, # fixed_τ
                200, # n_beads (if τ is not fixed)
                n_steps, # No. of steps
                50000, # number of thermalisation steps
                "Single", # movers
                "Frohlich", # potential type
                "Virial", # estimators
                false, # Not quick steps
                false, # threading
                1.0, # Start Range
                1,
                1, 
                0.005, # observable skips
                0.5, # equilibrium skips
                1,
                true,
                1
                );
end

begin
    # Define plot parameters
    default(fontfamily="Times New Roman",
        titlefont = (16, "Computer Modern"),
        guidefont = (18, "Computer Modern"),
        tickfont = (12, "Computer Modern"),
        legendfontsize = 12,
        linewidth=2, framestyle=:box, label=nothing, grid=true)

    # Energy/Position Plots
    positions_flatten = collect(Iterators.flatten(positions))
    posplot = histogram(positions_flatten, xlab = "Position")
    display(posplot)
    energy_plot = plot(energies, ylabel="Energy", xlab = "Sweeps")
    hline!([comparison_energy], linestyle=:dash)
    energy_hist = histogram(energies, ylab="Frequencies", xlab="Energy")
    positions_flatten = collect(Iterators.flatten(positions))
    #posplot = histogram(positions_flatten, xlab = "Position")

    # Correlation cut-off set by n
    #corr_plot = plot(1:n_beads-1, corr_mean[1:n_beads-1], yerror = corr_std, ylabel="G(Δτ)", xlabel = "Δτ")

    # Autocorrelation plot
    autoCorrelation1 = autoCorrelation(energies)
    auto_plot = plot(1:length(energies)-1, autoCorrelation1[1:length(energies)-1], ylabel=L"C_{k}", xlab = "k / $observables_skip\$ n\$")
    println("correlation time is:", autoCorrelationTime(autoCorrelation1))

    # Displaying plot command
    #display(energy_hist)
    display(energy_plot)

    #display(auto_plot)

    SaveJLDData(T, potential, version, n_beads, n_steps, data)
end
