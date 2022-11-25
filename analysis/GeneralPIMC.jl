using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility
using LaTeXStrings


function generalPIMC(T, m, ω, α, n_particles, n_dimensions, regime, fixing_tau, fixed_τ, n_beads, n_steps, n_thermalised, movers, potential, estimator, threads::Bool = true, start_range = 1.0)
    
    """
    Initialise System Variables
    """

    # Path variables
    β = 1 / T
    ħ = 1.0
            
    # For fixed τ 
    if fixing_tau
        adjusted_beads = Int(floor(1/(fixed_τ*T)))
        if adjusted_beads < 1
            adjusted_beads = 1
        end
        path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)

    # For fixed number of beads
    else
        τ = 1.0 / (T * n_beads)
        path = Path(n_beads, n_particles, n_dimensions, τ, m = m)
    end

    if threads
        plot_on = false
    else
        plot_on = true
    end
    """
    Set Potential Function
    """
    
    if potential == "Frohlich"
        potential = FrohlichPotential(α,ω,ħ) # ω is phonon frequency
    elseif potential == "Harmonic"
        potential = HarmonicPotential(ω)
    elseif potential == "MexicanHat"
        potential = MexicanHatPotential(80000.0)
    elseif potential == "Constant"
        potential = ConstantPotential(10.0)
    end

    """
    PIMC Variables
    """

    #skipping between sampling
    equilibrium_skip = 0.1 * n_steps
    observables_skip = 0.0002 * n_steps

    observables = [Energy, Position]
    
    if estimator == "Virial"
        estimators = [Virial_Estimator()]

    elseif estimator == "Thermodynamic"
        estimators = [Thermodynamic_Estimator()]
    
    elseif estimator == "Simple"
        estimators = [Simple_Estimator()]
    end
    """
    Run Simulation
    """

    thermalised_start!(path, potential, n_steps = n_thermalised)
    pimc = PIMC(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true, visual=false)
    acceptance_ratio = pimc[1]
    output_observables = pimc[2]

    energy = output_observables["Energy"][string(Symbol(estimators[1]))]
    position = output_observables["Position"][string(Symbol(estimators[1]))]

    # Comparison energy
    if typeof(potential) == HarmonicPotential
        comparison_energy = analytic_energy_harmonic(potential,β,ħ)
    elseif typeof(potential) == FrohlichPotential
        #comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true) # orginally threads is true
        #comparison_energy = comparison_polaron.F
        comparison_energy = 0.0
    end

    variances = jackknife(energy)

    # Post analysis
    
    println("acceptance ratio = ", acceptance_ratio)
    println("Mean energy = ", mean(energy))
    println("comparison_energy = ", comparison_energy)
    println("jackknife errors = ", sqrt(variances[2]))
    println("number of energy ", length(energy))
    

    #Plots
    
    if plot_on
        energyplot = scatter(energy[Int(floor(0.8*length(energy))):end], ylabel="Energy", xlabel="x * $observables_skip steps", labels = estimator)
        #energyplot = scatter(energy, ylabel="Energy", xlabel="x * $observables_skip steps")
        title!("Energy against steps")
        display(energyplot)
        savefig("./figs/Energy_conv/Energyvssteps_T=$(T)_energy=$(estimator).png")
    end
    
    #posplot = histogram(position[:,1,1])
    #plot(posplot, energyplot, layout = (2,1), legend = false)
    #plot(posplot, xlabel="Position", ylabel="Prob Amplitude", legend = false)

    return energy, variances, acceptance_ratio

end

#=
begin
    generalPIMC(0.1, #Temperature
                1.0, # mass
                1.0, # ω (has to be float)
                1.0, # α (has to be float)
                1, # no of particles
                1, # number of n_dimensions
                Primitive_Regime(), # regime type
                true, # fixing tau or not
                0.1, # fixed_τ
                200, # n_beads of tau not fixed
                100000, # No. of steps
                100000, # number of thermalisation
                [[Single!],[1.0]], # movers
                "Frohlich", # potential type
                "Virial", # estimators
                false # not threading
                )
end
=#
