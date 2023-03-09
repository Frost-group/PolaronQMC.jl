using Revise
using PolaronQMC
using Statistics
using Plots
using Dates
using PolaronMobility
using LaTeXStrings
using CSV

@time begin
    
    """
    Initialise System Variables
    """

    # Randomly select a version number
    version = Int(floor(rand()*1000))

    # Set Parameters
    T = 0.1
    m = 1.0
    ω = 1.0
    α = 1.0
    ħ = 1.0

    n_particles = 2
    n_dimensions = 3
    start_range = 1.0
    β = 1 / T

    # Number of Monte-Carlo-Steps
    n_steps = 500

    # Choose potential from "Harmonic", "Frohlich", "MexicanHat", "Constant"
    potential = "Harmonic"

    # Choose Monte-Carlo Mover from "Single", "Displace", "Bisect"
    mover = "Single"

    # Choose path regime from "Simple", "Primitive"
    regime = "Primitive"
    
    # Choose energy estimator from "Simple", "Virial", "Thermodynamic"
    estimator = "Virial"
    
    # Choose Observables
    energies = true
    positions = true
    
    # Choose Estimators
    energy_estimators = []

    # Pick True for fixed beads or False for fixed τ
    fixed_beads = true
    if fixed_beads
        n_beads = 1500
        τ = 1.0 / (T * n_beads)
    else
        # For fixed τ
        τ = 0.1
        n_beads = Int(floor(1/(fixed_τ*T)))
    end

    # Fixed observable skip or step dependant
    if true
        equilibrium_skip = 1 
        observables_skip = 1
    else
        equilibrium_skip = 0.5 * n_steps
        observables_skip = 0.001 * n_steps
    end

    # Initate path
    path = Path(n_beads, n_particles, n_dimensions, τ, m=m)

    # Set regime
    if regime == "Primitive"
        regime = PrimitiveRegime()
    elseif regime == "Simple"
        regime = SimpleRegime()
    else
        println("Invalid Regime: ", regime)
    end
   
    # Set Potential
    if potential == "Frohlich"
        potential = FrohlichPotential(α,ω,ħ)
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

    # thermalised_start!(path, potential, n_steps = 100000)
    data = PIMC(n_steps, equilibrium_skip, observables_skip, path, mover, estimators, potential, regime, energies, positions, adjust=true)
    
    # Store outputs
    energies = data["Energy:$(estimator)"]
    positions = data["Position:p1d1"]

    # Flatten position matrix to Array
    positions = collect(Iterators.flatten(positions))

    # Comparison energy
    if typeof(potential) == HarmonicPotential
        comparison_energy = analyticEnergyHarmonic(pot,β,ħ,n_dimensions)
    elseif typeof(potential) == FrohlichPotential
        comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
        comparison_energy = comparison_polaron.F
    end

    # Post analysis
    variances = jackknife(energies)
    jacknife_errors = sqrt(variances[2])
    mean_energy = mean(energies)

    # Output measurements and statistics
    println("Number of Beads: ", n_beads)
    println("α: ", α)
    println("Mean Energy: ", mean_energy)
    println("Comparison Energy: ", comparison_energy)
    println("jackknife errors: ", jacknife_errors)


    # Plots
    energy_plot = plot(energies, ylabel="Energy", xlab = "Sweeps / $observables_skip\$ n\$")
    posplot = histogram(positions)
    display(energy_plot)
    display(posplot)

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


