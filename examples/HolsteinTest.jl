# Holstein test

using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility
using LaTeXStrings
using DelimitedFiles
using JLD

@time begin
    
    """
    Initialise System Variables
    """

    # Randomly select a version number for saving purposes
    version = Int(floor(rand()*10000))

    # Set Parameters, all use atomic units where m = ħ = ω = 1.0
    T = 0.01
    m = 1.0
    ω = 1.0
    α = 2.2
    ħ = 1.0
    J = 1.0 #hopping integral

    n_particles = 1
    n_dimensions = 3
    start_range = 1.0
    β = 1 / T

    # Number of Monte-Carlo-Steps
    n_steps = 1000000

    # Choose potential from "Harmonic", "Frohlich", "MexicanHat", "Constant"
    potential = "Holstein"
    pot = potential
    
    # Choose energy estimator from "Simple", "Virial", "Thermodynamic"
    estimator = "Virial"
    
    # Choose Observables
    observables = ["Energy", "Position"]
    
    # Choose Estimators
    energy_estimators = []

    # Pick True for fixed beads or False for fixed τ
    fixed_beads = false
    if fixed_beads
        n_beads = 100
        τ = 1.0 / (T * n_beads)
    else
        # For fixed τ
        τ = 0.05
        n_beads = Int(floor(1/(τ*T)))
    end

    # Fixed observable skip or step dependant
    quick_steps = true
    if quick_steps
        equilibrium_skip = 1
        observables_skip = 100
    else
        equilibrium_skip = 0.5 * n_steps #try to put as 0.5, for 0.2 for quicker testing process
        observables_skip = 0.005 * n_steps
    end

    # Initate path
    path = DiscretePath(n_beads, n_particles, n_dimensions, τ, m)
   
    # Set Potential
    if potential == "Holstein"
        potential = HolsteinPotential(α, ω, ħ, J)
    end

    # Set Estimator
    if estimator == "Virial"
        estimators = [VirialEstimator()]
    elseif estimator == "Thermodynamic"
        estimators = [ThermodynamicEstimator()]
    elseif estimator == "Simple"
        estimators = [SimpleEstimator()]
    else
        println("Invalid Estimator: ", estimator)
    end

    """
    Run Simulation
    """

    println("started T is ", T)
    println("n_step is ", n_steps)
    println("---started----!")

    data = Holstein_PIMC(n_steps, equilibrium_skip, observables_skip, path, estimators, potential, observables)
    
    # Store outputs
    energies = data["Energy:$(estimator)"]
    positions = data["Position:p1d1"] # Select a particular particle and dimension
    #correlations = data["Correlation:$(estimator)"]

    # Flatten position matrix to Array
    positions_flatten = collect(Iterators.flatten(positions))

    # Post analysis
    variances = jackknife(energies)
    jacknife_errors = sqrt(variances[2])
    mean_energy = mean(energies)

    # Saving data in a big jld file (dictionary)
    #save("data_arr/$(pot)/$(string(Symbol(potential)))_T$(T)_nsteps$(n_steps)_v$(version)_beads$(n_beads).jld", "data", data, "energies", energies, "comparison_energy", comparison_energy, "correlations", correlations, "jacknife_errors", jacknife_errors, 
        #        "equilibrium_skip", equilibrium_skip, "observables_skip", observables_skip, "final_pos", path.beads[:, :, :])
    
    # Output measurements and statistics
    println("version is:", version)
    println("Number of Beads: ", n_beads)
    println("Number of steps: ", n_steps)
    println("τ is: ", τ)
    println("Temperature: ", T)
    println("α: ", α)
    println("Mean Energy: ", mean_energy)
    #println("Comparison Energy: ", comparison_energy)
    println("jackknife errors: ", jacknife_errors)


        # Plots
    # Define plot parameters
    default(fontfamily="Times New Roman",
        titlefont = (16, "Computer Modern"),
        guidefont = (18, "Computer Modern"),
        tickfont = (12, "Computer Modern"),
        legendfontsize = 12,
        linewidth=2, framestyle=:box, label=nothing, grid=true)

    posplot = histogram(positions_flatten, xlab = "Position")
    display(posplot)

    energy_plot = plot(energies, ylabel="Energy", xlab = "Sweeps / $observables_skip\$ n\$")
    #hline!([comparison_energy], linestyle=:dash)
    energy_hist = histogram(energies, ylab="Frequencies", xlab="Energy")
    
    # Displaying plot command
    display(energy_hist)
    display(energy_plot)

    # Visualise 
    #animation = animate_PIMC(data, 1)
    #gif(animation, "image/anim_output_$(λ)_version$(version)_1.gif", fps = 50) 


end


