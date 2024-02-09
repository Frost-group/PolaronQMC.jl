# GeneralTest.jl -> Test for Harmonic Oscillator

using Revise
using PolaronQMC
using Statistics
using Plots
using Dates
using PolaronMobility
using LaTeXStrings
using JLD

@time begin

    """
    Initialise System Variables
    """

    # Randomly select a version number for saving purposes
    version = Int(floor(rand() * 10000))

    # Set PIMC ring parameters
    T = 0.1
    β = 1 / T
    m = 1.0
    n_particles = 1
    n_dimensions = 1
    start_range = 1.0

    # Set potential parameters
    ω = 1.0
    α = 1.0
    ħ = 1.0

    # Number of Monte-Carlo-Steps
    n_steps = 500000

    # Choose potential from "Harmonic", "Frohlich", "MexicanHat", "Constant"
    potential = "Harmonic"
    pot = potential

    # Choose Monte-Carlo Mover from "Single", "Displace", "Bisect"
    mover = "Single"

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
        n_beads = 500
        τ = 1.0 / (T * n_beads)
    else
        τ = 0.5
        n_beads = Int(floor(1 / (τ * T)))
    end

    # Fixed observable skip or step dependant
    quick_steps = false
    if quick_steps
        equilibrium_skip = 1000
        observables_skip = 1000
    else
        equilibrium_skip = 0.2 * n_steps #try to put as 0.5, 0.2 is for quicker testing process. Convergence has to be checked.
        observables_skip = 0.001 * n_steps
    end

    # Initate path
    path = Path(n_beads, n_particles, n_dimensions, τ, m = m)

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
    elseif potential == "Constant"
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
        comparison_energy =
            analyticEnergyHarmonic(potential.ω, β, ħ, n_dimensions) * n_particles
    elseif typeof(potential) == FrohlichPotential
        comparison_polaron =
            make_polaron(
                [α],
                [T],
                [0.0];
                ω = 1.0,
                rtol = 1e-4,
                verbose = true,
                threads = true,
            ) * n_particles
        comparison_energy = comparison_polaron.F
    end

    # Post analysis
    variances = jackknife(energies)
    jacknife_errors = sqrt(variances[2])
    mean_energy = mean(energies)
    corr_mean = mean(correlations)
    corr_std = std(correlations)

    #last_acceptance_rate = last(acceptance_rates)
    #mean_acceptance_rate = mean(acceptance_rates)
    #std_acceptance_rate = std(acceptance_rates)

    save(
        "data_arr/Harmonic/$(string(Symbol(potential)))_T$(T)_nsteps$(n_steps)_v$(version)_beads$(n_beads).jld",
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
    )

    # Output measurements and statistics
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
        fontfamily = "Computer Modern",
        titlefont = (16, "Computer Modern"),
        guidefont = (18, "Computer Modern"),
        tickfont = (12, "Computer Modern"),
        legendfontsize = 12,
        linewidth = 2,
        framestyle = :box,
        label = nothing,
        grid = true,
    )

    # Plots (e.g. Energy, position, Correlation)
    energy_plot =
        plot(energies, ylabel = "Energy", xlab = "Sweeps / $observables_skip\$ n\$")
    hline!([comparison_energy], linestyle = :dash)
    energy_hist = histogram(energies, ylab = "Frequencies", xlab = "Energy")
    posplot = histogram(positions_flatten, xlab = "Position")

    n = Int(floor(n_beads - 1))
    corr_plot =
        plot(1:n, corr_mean[1:n], yerror = corr_std[1:n], ylabel = "G(Δτ)", xlabel = "Δτ")

    #acceptance_rate_plot = plot(acceptance_rates[Int(length(acceptance_rates)*0.9):end], xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Acceptance\, Rate\, /\, } r", dpi=600)
    #shift_width_plot = plot(shift_widths, xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Shift\, Width\, /\, } \Delta x", dpi=600)

    display(energy_hist)
    display(energy_plot)
    display(posplot)
    display(corr_plot)

    if n_particles == 2
        positions1 = positions
        positions2 = collect(Iterators.flatten(data["Position:p2d1"]))
        posplot = histogram([positions1, positions2])
        display(posplot)
    end

    # Visualise 
    #anim = animate_PIMC(pimc, n_particles, n_dimensions, "3D Harmonic Potential", "Single 1.0 Mover", "0.1")
    #gif(anim, "saved_plots/anim_output.gif", fps = 60) 

end
