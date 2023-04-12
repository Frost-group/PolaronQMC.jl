using DelimitedFiles
using JLD

begin
    #Setting up folder naming
    particle_index = 1
    dimension_index = 1
    estimator = "Virial"

    file = "data_arr/Frohlich/FrohlichPotential(5.0, 1.0, 1.0)_α5.0_T0.08_nsteps150000_v413_beads2500_τ0.005.jld"
    full_data = load(file)
    data = load(file)["data"] # Extracting the full data dictionary
    position1 = data["Position:p$(particle_index)d$(dimension_index)"]
    positions_flatten = collect(Iterators.flatten(position1))
    energies = data["Energy:$(estimator)"]
    #acceptance_rates = load(file)["acceptance_rates"]
    #shift_widths = load(file)["shift_widths"]
    jacknife_errors = full_data["jacknife_errors"]
    comparison_energy = full_data["comparison_energy"]

    # Plots
    energy_plot = plot(energies, ylabel="Energy", xlab = "Sweeps")
    hline!([comparison_energy], linestyle=:dash)
    energy_hist = histogram(energies, ylab="Frequencies", xlab="Energy")
    posplot = histogram(positions_flatten, xlab = "Position")

    #=
    acceptance_rate_plot = plot(acceptance_rates[Int(length(acceptance_rates)*0.9):end], xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Acceptance\, Rate\, /\, } r", dpi=600)
    shift_width_plot = plot(shift_widths, xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Shift\, Width\, /\, } \Delta x", dpi=600)
    acceptance_shift_plot = scatter(acceptance_rates, shift_widths, xlab=L"\mathrm{Acceptance\, Rate\, /\, } r", ylab=L"\mathrm{Shift\, Width\, /\, } \Delta x", dpi=600)
    acceptance_rate_hist = histogram(acceptance_rates, ylab="Frequency", xlab=L"\mathrm{Acceptance\, Rate\, /\, } r")
    shift_width_hist = histogram(shift_widths, ylab="Frequency", xlab=L"\mathrm{Shift\, Width\, /\, } \Delta x")
    =#

    # Command in VSCode to display graphs
    display(energy_hist)
    display(energy_plot)
    display(posplot)
    #display(acceptance_rate_plot)
    #display(shift_width_plot)
    #display(acceptance_shift_plot)
    #display(shift_width_hist)
    #display(acceptance_rate_hist)


    println("Mean Energy: ", mean(energies))
    println("Comparison Energy: ", comparison_energy)
    println("jackknife errors: ", jacknife_errors)
    #println("Final Acceptance Rate: ", last_acceptance_rate)
    #println("Mean Acceptance Rate: ", mean_acceptance_rate, " +/- ", std_acceptance_rate)

end