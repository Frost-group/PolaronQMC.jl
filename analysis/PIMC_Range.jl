using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility
using LaTeXStrings
using DelimitedFiles
using Base.Threads
using JLD

function RangeAlphaPIMC(T::Float64, m::Float64, ω::Float64, α_array::Vector{Float64}, n_particles::Int64, n_dimensions::Int64, regime::String, fixed_beads::Bool, fixed_τ::Float64, n_beads::Int64, 
    n_steps::Int64, n_thermalised::Int64, mover::String, potential::String, estimator::String; quick_steps=false, threads::Bool = false, start_range = 1.0, particleIndex = 1, dimensionIndex = 1,
    observable_skip_factor=0.005, equilibrium_skip_factor=0.5, version = 1, verbose::Bool = true, thread_number::Int64 = 16)
    
    """
    Allow array of Alpha values.

    Initialise System Variables:
    T::Float or Array
    m::Float
    n_particles::Int
    n_dimensions::Int
    start_range::Float

    Initialise potential variables:
    ω::Float
    α::Float or Array of values

    Possible choice of potential: {"Harmonic", "Frohlich", "MexicanHat", "Constant"}
    Possible choice of estimator: {"Simple", "Virial", "Thermodynamic"}
    Possible choice of regime: {"Simple", "Primitive"}
    Possible choice of mover: {"Single", "Displace", "Bisect"}
    """
    
    # Storing energy values & Jackknife errors for different alpha values
    aggregated_data = Dict() # Create an empty dictionary for storing data
    Mean_energy_arr = [0.0 for j in 1:length(α_array)]
    Error_arr = [0.0 for j in 1:length(α_array)]

    β = 1 / T
    ħ = 1.0

    if fixed_beads
        n_beads = n_beads
        τ = 1.0 / (T * n_beads)
    else
        τ = fixed_τ
        n_beads = Int(floor(1/(τ*T)))
    end

    #Initialsing path
    path = Path(n_beads, n_particles, n_dimensions, τ, m=m)

    if threads
        plot_on = false
    else
        plot_on = true
    end

    for i in 1:length(α_array)

        println("Simulating multiple alphas")
        α = α_array[i]

        # Set potential function
        if potential == "Frohlich"
            potential = FrohlichPotential(α,ω,ħ) # ω is phonon frequency
        elseif potential == "Harmonic"
            potential = HarmonicPotential(ω)
        elseif potential == "MexicanHat"
            potential = MexicanHatPotential(80000.0)
        elseif potential == "Constant"
            potential = ConstantPotential(10.0)
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

        # Set Regime
        if regime == "Primitive"
            regime = PrimitiveRegime()
        elseif regime == "Simple"
            regime = SimpleRegime()
        elseif regime == "LBRegime"
            regime = LBRegime()
        else
            println("Invalid Regime: ", regime)
        end

        # Observables type
        observables = ["Energy", "Position", "Correlation"]

        # Fixed observable skip or step dependant
        if quick_steps
            equilibrium_skip = 10
            observables_skip = 1000
        else
            equilibrium_skip = equilibrium_skip_factor * n_steps #try to put as 0.5, 0.2 is for quicker testing process. Convergence has to be checked.
            observables_skip = observable_skip_factor * n_steps # How many steps to skip before making observables
        end

        """
        Run Simulation
        """

        println("------started--------")
        println("Temperature is: ", T)
        println("α is: ", α)
        println("n_step is: ", n_steps)

        # Thermalisation Process (With simple Potential like Harmonic)
        PIMC(n_thermalised, n_thermalised, n_thermalised, path, SingleMover(path), estimators, HarmonicPotential(1.0), regime, observables, adjust=true)
        if verbose
            println("Thermalisation complete")
        end

        #data = PIMC(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true, visual=false)
        data = PIMC(n_steps, equilibrium_skip, observables_skip, path, mover, estimators, potential, regime, observables, adjust=true)
        
        # Storing all the different outputz
        energies = data["Energy:$(estimator)"]
        positions = data["Position:p$(particleIndex)d$(dimensionIndex)"] # Select a particular particle and dimension
        correlations = data["Correlation:$(estimator)"]
        acceptance_rates = data["Acceptance Rate:p$(particleIndex)"]
        adjuster_values = data["Adjuster Value:p$(particleIndex)"]

        # Comparison energy
        if typeof(potential) == HarmonicPotential
            comparison_energy = analyticEnergyHarmonic(potential.ω,β,ħ,n_dimensions)
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
        println("Number of steps: ", n_steps)
        println("τ is: ", τ)
        println("Temperature: ", T)
        println("α: ", α)
        println("Mean Energy: ", mean_energy)
        println("Comparison Energy: ", comparison_energy)
        println("jackknife errors: ", jacknife_errors)
        
        Mean_energy_arr[i] = mean_energy
        Error_arr[i] = jacknife_errors
        aggregated_data[i] = data
    end

    return mean_energy, jacknife_errors, comparison_energy, energies, positions, correlations, acceptance_rates, adjuster_values, equilibrium_skip, observables_skip, version, potential, n_beads, data
end

function RangeTempPIMC(T_array::Vector{Float64}, m::Float64, ω::Float64, α::Float64, n_particles::Int64, n_dimensions::Int64, regime::String, fixed_beads::Bool, fixed_τ::Float64, n_beads::Int64, 
    n_steps::Int64, n_thermalised::Int64, mover::String, potential::String, estimator::String; quick_steps=false, threads::Bool = false, start_range = 1.0, particleIndex = 1, dimensionIndex = 1,
    observable_skip_factor=0.005, equilibrium_skip_factor=0.5, version = 1, verbose::Bool = true, thread_number::Int64 = 16)
    
    """
    Allow array of Temperature values.

    Initialise System Variables:
    T::Array of values
    m::Float
    n_particles::Int
    n_dimensions::Int
    start_range::Float

    Initialise potential variables:
    ω::Float
    α::Float

    Possible choice of potential: {"Harmonic", "Frohlich", "MexicanHat", "Constant"}
    Possible choice of estimator: {"Simple", "Virial", "Thermodynamic"}
    Possible choice of regime: {"Simple", "Primitive"}
    Possible choice of mover: {"Single", "Displace", "Bisect"}
    """

    # Storing energy values & Jackknife errors for different Temperature
    # Temperature changes the number of beads
    aggregated_data = Dict() # Create an empty dictionary for storing data
    Mean_energy_arr = [0.0 for j in 1:length(T_array)]
    Error_arr = [0.0 for j in 1:length(T_array)]

    # Defining parameters unaffected by temperature
    ħ = 1.0

    # Set potential function
    if potential == "Frohlich"
        potential = FrohlichPotential(α,ω,ħ) # ω is phonon frequency
    elseif potential == "Harmonic"
        potential = HarmonicPotential(ω)
    elseif potential == "MexicanHat"
        potential = MexicanHatPotential(80000.0)
    elseif potential == "Constant"
        potential = ConstantPotential(10.0)
    end

    if threads
        plot_on = false
    else
        plot_on = true
    end

    println("Simulating multiple temperatures")
    # Big for-loop, multi-threading possible
    @threads for i in 1:length(T_array)
        println("i = $i on thread $(Threads.threadid())")
        T = T_array[i]
        β = 1 / T

        if fixed_beads
            n_beads = n_beads
            τ = 1.0 / (T * n_beads)
        else
            τ = fixed_τ
            n_beads = Int(floor(1/(τ*T)))
        end

        #Initialsing path (Temperature-dependent)
        path = Path(n_beads, n_particles, n_dimensions, τ, m=m)

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

        # Set Regime
        if regime == "Primitive"
            regime = PrimitiveRegime()
        elseif regime == "Simple"
            regime = SimpleRegime()
        elseif regime == "LBRegime"
            regime = LBRegime()
        else
            println("Invalid Regime: ", regime)
        end

        # Observables type
        observables = ["Energy", "Position", "Correlation"]

        # Fixed observable skip or step dependant
        if quick_steps
            equilibrium_skip = 10
            observables_skip = 1000
        else
            equilibrium_skip = equilibrium_skip_factor * n_steps #try to put as 0.5, 0.2 is for quicker testing process. Convergence has to be checked.
            observables_skip = observable_skip_factor * n_steps # How many steps to skip before making observables
        end

        """
        Run Simulation
        """

        println("------started--------")
        println("Temperature is: ", T)
        println("α is: ", α)
        println("n_step is: ", n_steps)

        # Thermalisation Process (With simple Potential like Harmonic)
        PIMC(n_thermalised, n_thermalised, n_thermalised, path, SingleMover(path), estimators, HarmonicPotential(1.0), regime, observables, adjust=true)
        if verbose
            println("Thermalisation complete")
        end

        #data = PIMC(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true, visual=false)
        data = PIMC(n_steps, equilibrium_skip, observables_skip, path, mover, estimators, potential, regime, observables, adjust=true)
        
        # Comparison energy
        if typeof(potential) == HarmonicPotential
            comparison_energy = analyticEnergyHarmonic(potential.ω,β,ħ,n_dimensions)
        elseif typeof(potential) == FrohlichPotential
            comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
            comparison_energy = comparison_polaron.F
        end

        # Storing all the different outputz
        energies = data["Energy:$(estimator)"]
        positions = data["Position:p$(particleIndex)d$(dimensionIndex)"] # Select a particular particle and dimension
        correlations = data["Correlation:$(estimator)"]
        acceptance_rates = data["Acceptance Rate:p$(particleIndex)"]
        adjuster_values = data["Adjuster Value:p$(particleIndex)"]
        data["Comparsion Energy"] = comparison_energy


        # Post analysis
        variances = jackknife(energies)
        jacknife_errors = sqrt(variances[2])
        mean_energy = mean(energies)

        # Output measurements and statistics
        println("Number of Beads: ", n_beads)
        println("Number of steps: ", n_steps)
        println("τ is: ", τ)
        println("Temperature: ", T)
        println("α: ", α)
        println("Mean Energy: ", mean_energy)
        println("Comparison Energy: ", comparison_energy)
        println("jackknife errors: ", jacknife_errors)
        
        Mean_energy_arr[i] = mean_energy
        Error_arr[i] = jacknife_errors
        aggregated_data[i] = data
    end

    return aggregated_data, Mean_energy_arr, Error_arr
end

function SaveJLDData(T, potential, version, n_beads, n_steps, data)
    save("data_arr/$(string(typeof(potential)))/$(string(Symbol(potential)))_T$(T)_nsteps$(n_steps)_v$(version)_beads$(n_beads).jld", "data", data)
end