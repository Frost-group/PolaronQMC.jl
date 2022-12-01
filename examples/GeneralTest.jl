using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility
using LaTeXStrings


begin
    
    """
    Initialise System Variables
    """

    # Path variables
    T = 0.1
    m = 1.0
    n_particles = 1
    n_dimensions = 3
    start_range = 1.0
    β = 1 / T

    # For fixed τ 
    fixed_τ = 0.1
    adjusted_beads = Int(floor(1/(fixed_τ*T)))

    # For fixed number of beads
    n_beads = 200
    τ = 1.0 / (T * n_beads)

    # path = Path(n_beads, n_particles, n_dimensions, τ, m = m)
    path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)

    # Set regime
    regime = Primitive_Regime()

    """
    Set Potential Function
    """
    
    # Potential variables
    ω = 1.0
    α = 1.0
    ħ = 1.0
    
    #potential = FrohlichPotential(α,ω,ħ)
    potential = HarmonicPotential(ω)
    #potential = MexicanHatPotential(80.0)
    #potential = ConstantPotential(10.0)

    """
    PIMC Variables
    """

    #number of steps
    n_steps = 10000

    #skipping between sampling
    #equilibrium_skip = 0.5 * n_steps
    equilibrium_skip = 0
    observables_skip = 0.001 * n_steps
    #observables_skip = 1

    #types of moves
    #movers = [[Bisect!],[1.0]]
    #movers = Dict("Single!" => [1.0])
    movers = Dict("Displace" => [1.0])
    #movers = Dict("Single!" => [1.0], "Displace!" => [0.2])


    observables = [Energy, Position]
    
    estimators = [Virial_Estimator()]
    #estimators = [Thermodynamic_Estimator()]
    #estimators = [Simple_Estimator()]
 
    """
    Run Simulation
    """

    # thermalised_start!(path, potential, n_steps = 100000)
    pimc = PIMC(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true, visual=true)
    
    # Store outputs
    adjuster_stats = pimc[1]
    output_observables = pimc[2]

    energies = output_observables["Energy"][string(Symbol(estimators[1]))]
    acceptance_rates = adjuster_stats["Single!"]["Acceptance Rate"]
    shift_widths = adjuster_stats["Single!"]["Shift Width"]
    # position = output_observables["Position"][string(Symbol(estimators[1]))]

    # Comparison energy
    if typeof(potential) == HarmonicPotential
        comparison_energy = analytic_energy_harmonic(potential,β,ħ)
    elseif typeof(potential) == FrohlichPotential
        comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
        comparison_energy = comparison_polaron.F
    end

    # Post analysis
    variances = jackknife(energies)
    jacknife_errors = sqrt(variances[2])
    mean_energy = mean(energies)
    last_acceptance_rate = last(acceptance_rates)
    mean_acceptance_rate = mean(acceptance_rates)
    std_acceptance_rate = std(acceptance_rates)

    # Output measurements and statistics
    println("Number of Beads: ", adjusted_beads)
    println("Mean energy: ", mean_energy)
    println("Comparison Energy: ", comparison_energy)
    println("jackknife errors: ", jacknife_errors)
    println("Final Acceptance Rate: ", last_acceptance_rate)
    println("Mean Acceptance Rate: ", mean_acceptance_rate, " +/- ", std_acceptance_rate)

    # Define plot parameters
    default(fontfamily="Computer Modern",
        titlefont = (16, "Computer Modern"),
        guidefont = (18, "Computer Modern"),
        tickfont = (12, "Computer Modern"),
        legendfontsize = 12,
        linewidth=2, framestyle=:box, label=nothing, grid=true)

    # Plots
    energyplot = plot(energies, ylabel="Energy", xlab = "Sweeps / $observables_skip\$ n\$")
    acceptance_rate_plot = plot(acceptance_rates, xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Acceptance\, Rate\, /\, } r", dpi=600)
    shift_width_plot = plot(shift_widths, xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Shift\, Width\, /\, } \Delta x", dpi=600)
    acceptance_shift_plot = scatter(acceptance_rates, shift_widths, xlab=L"\mathrm{Acceptance\, Rate\, /\, } r", ylab=L"\mathrm{Shift\, Width\, /\, } \Delta x", dpi=600)

    #posplot = histogram(position[:,1,1])

    #plot(posplot, energyplot, layout = (2,1), legend = false)
    #plot(posplot, xlabel="Position", ylabel="Prob Amplitude", legend = false)

    display(energyplot)
    display(acceptance_rate_plot)
    display(shift_width_plot)
    display(acceptance_shift_plot)


    #savefig(acceptance_rate_plot, "saved_plots/acceptance_rate_convergence.png")

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

