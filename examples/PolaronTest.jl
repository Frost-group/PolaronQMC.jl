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
    version = Int(floor(rand() * 10000))

    # Set Parameters, all use atomic units where m = ħ = ω = 1.0
    T = 0.1
    m = 1.0
    ω = 1.0
    α = 3.0
    ħ = 1.0

    n_particles = 1
    n_dimensions = 3
    start_range = 1.0
    β = 1 / T

    # Number of Monte-Carlo-Steps
    n_steps = 20000

    # Choose potential from "Harmonic", "Frohlich", "MexicanHat", "Constant"
    potential = "Frohlich"
    pot = potential

    # Choose Monte-Carlo Mover from "Single", "Displace", "Bisect"
    mover = "Single"
    #mover = "Bisect"

    # Choose path regime from "Simple", "Primitive"
    regime = "Primitive"

    # Choose energy estimator from "Simple", "Virial", "Thermodynamic"
    estimator = "Virial"

    # Choose Observables
    observables = ["Energy", "Position", "Correlation"]

    # Choose Estimators
    energy_estimators = []

    # Pick True for fixed beads or False for fixed τ
    fixed_beads = false
    if fixed_beads
        n_beads = 250
        τ = 1.0 / (T * n_beads)
    else
        # For fixed τ
        τ = 0.01
        n_beads = Int(floor(1 / (τ * T)))
    end

    # Fixed observable skip or step dependant
    quick_steps = true
    if quick_steps
        equilibrium_skip = 10
        observables_skip = 100
    else
        equilibrium_skip = 0.5 * n_steps #try to put as 0.5, for 0.2 for quicker testing process
        observables_skip = 0.002 * n_steps
    end

    # Initate path
    path = Path(n_beads, n_particles, n_dimensions, τ, m = m, start_range = 10.0)

    # Set regime
    if regime == "Primitive"
        regime = PrimitiveRegime()
    elseif regime == "Simple"
        regime = SimpleRegime()
    elseif regime == "LBRegime"
        regime = LBRegime()
    else
        println("Invalid Regime: ", regime)
    end

    # Set Potential
    if potential == "Frohlich"
        potential = FrohlichPotential(α, ω, ħ)
    elseif potential == "Harmonic"
        potential = HarmonicPotential(ω)
    elseif potential == "MexicanHat"
        potential = MexicanHatPotential(80.0)
    elseif potential == "Contsant"
        potential = ConstantPotential(10.0)
    else
        println("Invalid Potential: ")
    end

    # Set Monte-Carlo Mover
    if mover == "Single"
        mover = SingleMover(path)
    elseif mover == "Displace"
        mover = DisplaceMover(path)
    elseif mover == "Bisect"
        mover = BisectMover(path)
    else
        println("Invalid mover")
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
    println("started α is ", α)
    println("n_step is ", n_steps)
    # thermalised_start!(path, potential, n_steps = 100000)
    data = PIMC(
        n_steps,
        equilibrium_skip,
        observables_skip,
        path,
        mover,
        estimators,
        potential,
        regime,
        observables,
        adjust = true,
    )

    # Store outputs
    energies = data["Energy:$(estimator)"]
    positions = data["Position:p1d1"] # Select a particular particle and dimension
    correlations = data["Correlation:$(estimator)"]

    # Flatten position matrix to Array
    positions_flatten = collect(Iterators.flatten(positions))

    # Comparison energy
    if typeof(potential) == HarmonicPotential
        comparison_energy = analyticEnergyHarmonic(potential, β, ħ, n_dimensions)
    elseif typeof(potential) == FrohlichPotential
        comparison_polaron = make_polaron(
            [α],
            [T],
            [0.0];
            ω = 1.0,
            rtol = 1e-4,
            verbose = true,
            threads = true,
        )
        comparison_energy = comparison_polaron.F
    end

    # Post analysis
    variances = jackknife(energies)
    jacknife_errors = sqrt(variances[2])
    mean_energy = mean(energies)
    corr_mean = mean(correlations)
    corr_std = std(correlations)

    # Saving data in a big jld file (dictionary)
    save(
        "data_arr/$(pot)/$(string(Symbol(potential)))_T$(T)_nsteps$(n_steps)_v$(version)_beads$(n_beads).jld",
        "data",
        data,
        "energies",
        energies,
        "comparison_energy",
        comparison_energy,
        "correlations",
        correlations,
        "jacknife_errors",
        jacknife_errors,
        "equilibrium_skip",
        equilibrium_skip,
        "observables_skip",
        observables_skip,
        "final_pos",
        path.beads[:, :, :],
    )


    # Output measurements and statistics
    println("version is:", version)
    println("Number of Beads: ", n_beads)
    println("Number of steps: ", n_steps)
    println("τ is: ", τ)
    println("Temperature: ", T)
    println("α: ", α)
    println("Mean Energy: ", mean_energy)
    println("Comparison Energy: ", comparison_energy)
    println("jackknife errors: ", jacknife_errors)
    #println("Final Acceptance Rate: ", last_acceptance_rate)
    #println("Mean Acceptance Rate: ", mean_acceptance_rate, " +/- ", std_acceptance_rate)

    # Define plot parameters
    default(
        fontfamily = "Times New Roman",
        titlefont = (16, "Computer Modern"),
        guidefont = (18, "Computer Modern"),
        tickfont = (12, "Computer Modern"),
        legendfontsize = 12,
        linewidth = 2,
        framestyle = :box,
        label = nothing,
        grid = true,
    )

    # Plots
    energy_plot =
        plot(energies, ylabel = "Energy", xlab = "Sweeps / $observables_skip\$ n\$")
    hline!([comparison_energy], linestyle = :dash)
    energy_hist = histogram(energies, ylab = "Frequencies", xlab = "Energy")
    posplot = histogram(positions_flatten, xlab = "Position")

    # Correlation cut-off set by n
    n = n_beads - 1
    corr_plot =
        plot(1:n, corr_mean[1:n], yerror = corr_std, ylabel = "G(Δτ)", xlabel = "Δτ")

    #acceptance_rate_plot = plot(acceptance_rates[Int(length(acceptance_rates)*0.9):end], xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Acceptance\, Rate\, /\, } r", dpi=600)
    #shift_width_plot = plot(shift_widths, xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Shift\, Width\, /\, } \Delta x", dpi=600)

    # Displaying plot command
    display(energy_hist)
    display(energy_plot)
    display(posplot)
    #display(corr_plot)

    if n_particles == 2
        positions1 = positions
        positions2 = collect(Iterators.flatten(data["Position:p2d1"]))
        posplot = histogram([positions1, positions2])
        display(posplot)
    end

    #savefig(energy_hist, "saved_plots/energy_hist.png") 
    #savefig(energy_plot, "saved_plots/energy_plot.png")    
    #savefig(acceptance_rate_plot, "saved_plots/acceptance_rate_convergence_02_12.png")
    #savefig(shift_width_plot, "saved_plots/shift_width_convergence_02_12.png")
    #savefig(acceptance_shift_plot, "saved_plots/acceptance_shift_02_12.png")
    #savefig(shift_width_hist, "saved_plots/shift_width_hist_02_12.png")
    #savefig(acceptance_rate_hist, "saved_plots/acceptance_rate_shift_02_12.png")

    #=
    acceptanceshiftplot = scatter(adjuster_stats["Single!"]["Acceptance Rate"],
                            adjuster_stats["Single!"]["Shift Width"],
                            xlabel= "Acceptance Rate",
                            ylabel= "Shift Width")
    display(acceptanceshiftplot) =#

    # Visualise 
    #anim = animate_PIMC(pimc, n_particles, n_dimensions, "3D Harmonic Potential", "Single 1.0 Mover", "0.1")
    #gif(anim, "saved_plots/anim_output.gif", fps = 60) 

end
