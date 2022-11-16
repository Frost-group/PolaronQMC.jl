using Revise
using PolaronQMC
using Statistics
using Plots
include("../src/PolaronQMCVisualisation.jl")
using .PolaronQMCVisualisation
using PolaronMobility
using LaTeXStrings


function generalPIMC(T, m, n_particles, n_dimensions, regime, fixing_tau, fixed_τ, n_beads, n_steps, n_thermalised, movers, start_range = 1.0)
    
    """
    Initialise System Variables
    """

    # Path variables
    β = 1 / T
            
    # For fixed τ 
    if fixing_tau
        adjusted_beads = Int(floor(1/(fixed_τ*T)))
        path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)

    # For fixed number of beads
    else
        τ = 1.0 / (T * n_beads)
        path = Path(n_beads, n_particles, n_dimensions, τ, m = m)
    end

    """
    Set Potential Function
    """
    
    # Potential variables
    ω = 1.0
    α = 1.0
    ħ = 1.0
    
    #potential = FrohlichPotential(α,ω,ħ)
    potential = HarmonicPotential(ω)
    #potential = MexicanHatPotential(80000.0)
    #potential = ConstantPotential(10.0)

    """
    PIMC Variables
    """

    #skipping between sampling
    # equilibrium_skip = 0.5 * n_steps
    equilibrium_skip = 0.1 * n_steps
    #equilibrium_skip = 0
    #observables_skip = 0.001 * n_steps
    observables_skip = 500
    #observables_skip = 10 * n_steps

    

    observables = [Energy, Position]
    
    estimators = [Virial_Estimator()]
    #estimators = [Thermodynamic_Estimator()]
    #estimators = [Simple_Estimator()]
        
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
        comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
        comparison_energy = comparison_polaron.F
    end

    variances = jackknife(energy)

    return energy, variances
end


#types of moves
#movers = [[Bisect!],[1.0]]
movers = [[Single!],[1.0]]
#movers = [[Displace!,Single!],[0.2,1.0]]
generalPIMC(0.1, #Temperature
            1, # mass
            1, # no of particles
            1, # number of n_dimensions
            Simple_Regime(), # regime type
            true, # fixing tau or not
            0.4, # fixed_τ
            200, #n_beads of tau not fixed
            10000000, # No. of steps
            100000, # number of thermalisation
            movers
            )
