begin
    using Revise
    using PolaronQMC
    using Statistics
    using Plots
    using PolaronMobility
    using LaTeXStrings
    using BenchmarkTools
end

@time begin
    n_steps = 500000;
    est = "Thermodynamic"
    #est = "Thermodynamic"
    pot = "Harmonic";
    T = 0.1; version=rand(1:10000)
    data, energy, error, comparison_energy, equilibrium_skip, observables_skip, n_beads, path = 
        generalPIMC(T, 1.0, 2.0, 3, 500000, version=1, pot=pot, estimator=est, fixed_τ=0.2)
    println("-------simulation end------------")
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
    positions = data["Position:p1d1"]
    energies = data["Energy:$(est)"]

    positions_flatten = collect(Iterators.flatten(positions))
    posplot = histogram(positions_flatten, xlab = "Position")
    display(posplot)
    energy_plot = plot(energies, ylabel="Energy", xlab = "Sweeps")
    hline!([comparison_energy], linestyle=:dash)
    energy_hist = histogram(energies, ylab="Frequencies", xlab="Energy")
    positions_flatten = collect(Iterators.flatten(positions))

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

    #SaveJLDData(T, potential, version, n_beads, n_steps, data)
end
