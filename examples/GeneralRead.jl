using DelimitedFiles
using JLD

begin
    #file = "data_arr/$(string(Symbol(potential)))_shift_widths_α$(α)_nsteps$(n_steps).jld"
    file = "data_arr/Frohlich/Frohlich_T=0.1_α4.0_nsteps25000_dim=3_v357.jld"
    position1 = load(file)["position"]
    energies = load(file)["energies"]
    acceptance_rates = load(file)["acceptance_rates"]
    shift_widths = load(file)["shift_widths"]
    jacknife_errors = load(file)["jacknife_errors"]
    comparison_energy = load(file)["comparison_energy"]

    # Plots
    energy_plot = plot(energies, ylabel="Energy", xlab = "Sweeps / $observables_skip\$ n\$")
    hline!([comparison_energy], linestyle=:dash)
    energy_hist = histogram(energies, ylab="Frequencies", xlab="Energy")
    acceptance_rate_plot = plot(acceptance_rates[Int(length(acceptance_rates)*0.9):end], xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Acceptance\, Rate\, /\, } r", dpi=600)
    shift_width_plot = plot(shift_widths, xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Shift\, Width\, /\, } \Delta x", dpi=600)
    acceptance_shift_plot = scatter(acceptance_rates, shift_widths, xlab=L"\mathrm{Acceptance\, Rate\, /\, } r", ylab=L"\mathrm{Shift\, Width\, /\, } \Delta x", dpi=600)
    acceptance_rate_hist = histogram(acceptance_rates, ylab="Frequency", xlab=L"\mathrm{Acceptance\, Rate\, /\, } r")
    #shift_width_hist = histogram(shift_widths, ylab="Frequency", xlab=L"\mathrm{Shift\, Width\, /\, } \Delta x")
    #posplot = histogram(position[:,1,1])
    #plot(posplot, energyplot, layout = (2,1), legend = false)
    #plot(posplot, xlabel="Position", ylabel="Prob Amplitude", legend = false)
    display(energy_hist)
    display(energy_plot)
    display(acceptance_rate_plot)
    display(shift_width_plot)
    display(acceptance_shift_plot)
    #display(shift_width_hist)
    display(acceptance_rate_hist)
    posplot = histogram(position1[:,1,1], xlab = "Position")
    display(posplot)

    println("Mean Energy: ", mean(energies))
    println("Comparison Energy: ", comparison_energy)
    println("jackknife errors: ", jacknife_errors)
    println("Final Acceptance Rate: ", last_acceptance_rate)
    println("Mean Acceptance Rate: ", mean_acceptance_rate, " +/- ", std_acceptance_rate)

end