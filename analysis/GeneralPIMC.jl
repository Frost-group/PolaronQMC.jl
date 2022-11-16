using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility


function generalPIMC(T, m, n_particles, n_dimensions, start_range)
    
    """
    Initialise System Variables
    """

    # Path variables
    T = 0.1
    m = 1.0
    n_particles = 1
    n_dimensions = 1
    start_range = 1.0
    β = 1 / T
            
    # For fixed τ 
    fixed_τ = 0.01
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
    #potential = MexicanHatPotential(80000.0)
    #potential = ConstantPotential(10.0)

    """
    PIMC Variables
    """

    #number of steps
    n_steps = 200000

    #skipping between sampling
    equilibrium_skip = 0.5 * n_steps
    #equilibrium_skip = 0
    observables_skip = 0.001 * n_steps
    #observables_skip = 10 * n_steps

    #types of moves
    #movers = [[Bisect!],[1.0]]
    movers = [[Single!],[1.0]]
    #movers = [[Displace!,Single!],[0.2,1.0]]

    observables = [Energy, Position]
    
    estimators = [Virial_Estimator()]
    #estimators = [Thermodynamic_Estimator()]
    #estimators = [Simple_Estimator()]
        
    """
    Run Simulation
    """

    #thermalised_start!(path, potential, n_steps = 100000)
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


generalPIMC()
