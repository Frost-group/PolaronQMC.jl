using DelimitedFiles
using JLD

function ReadingSavedJLDFiles(file_path, particleIndex, dimensionIndex, estimator, verbose)
    """
    Required info:
    1. file_path (in string)
    2. particleIndex (int; for plotting position)
    3. dimensionIndex (int; for plotting position)
    4. estimator: {"Virial", "Thermodynamic", "Simple"}
    5. verbose (bool; whether to print out)
    """

    # Getting the full data dictionary
    full_data = load(file)
    data = load(file)["data"] # Extracting the full data dictionary, but printing the data can view all entries/keys

    # Getting positions and flattened to array
    position1 = data["Position:p$(particleIndex)d$(dimensionIndex)"]
    positions_flatten = collect(Iterators.flatten(position1))

    # Getting energies and uncertainties
    energies = data["Energy:$(estimator)"]
    jacknife_errors = full_data["jacknife_errors"]
    comparison_energy = full_data["comparison_energy"]

    # Plotting (to uncomment)
    #=
    energy_plot = plot(energies, ylabel="Energy", xlab = "Sweeps")
    hline!([comparison_energy], linestyle=:dash)
    energy_hist = histogram(energies, ylab="Frequencies", xlab="Energy")
    posplot = histogram(positions_flatten, xlab = "Position")

    # Command in VSCode to display graphs
    display(energy_hist)
    display(energy_plot)
    display(posplot)
    =#

    if verbose
        println("Mean Energy: ", mean(energies))
        println("Comparison Energy: ", comparison_energy)
        println("jackknife errors: ", jacknife_errors)
    end

    return energies, jacknife_errors, comparison_energy
end


begin
    file = "data_arr/Harmonic/HarmonicPotential(1.0)_T0.1_nsteps100000_v9912_beads5.jld"
    E, j, C = ReadingSavedJLDFiles(file, 1, 1, "Virial", true)
end
