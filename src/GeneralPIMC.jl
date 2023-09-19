using Revise
using Statistics
using Plots
using PolaronMobility
using LaTeXStrings
using JLD
using Base.Threads

function generalPIMC(T::Float64, ω::Union{Float64, Vector{Float64}}, α::Float64, n_dimensions::Int64, n_steps::Int64; 
    m::Float64=1.0, n_particles::Int64=1, regime::String="Primitive", fixed_beads::Bool=false, 
    fixed_τ::Float64=0.02, n_beads::Int64=50, n_thermalised::Int64=10000, mover::String="Single", 
    pot::String="Frohlich", estimator::String="Virial", quick_steps::Bool=false, threads::Bool = false, 
    start_range::Float64 = 1.0, particleIndex::Int64 = 1, dimensionIndex::Int64 = 1,
    observable_skip_factor::Float64=0.005, equilibrium_skip_factor::Float64=0.5, version::Int64 = 1, 
    verbose::Bool = true, thread_number::Int64 = 16)
    
    """
    Initialise System Variables:
    T
    m
    n_particles
    n_dimensions
    start_range

    Initialise potential variables:
    ω
    α

    Possible choice of potential: {"Harmonic", "Frohlich", "MexicanHat", "Constant"}
    Possible choice of estimator: {"Simple", "Virial", "Thermodynamic"}
    Possible choice of regime: {"Simple", "Primitive"}
    Possible choice of mover: {"Single", "Displace", "Bisect"}
    """
    
    #version = Int(floor(rand()*10000))

    # Path variables
    β = 1 / T
    ħ = 1.0

    # For fixed number of beads or not fixing beads
    if fixed_beads
        n_beads = n_beads
        #τ = 1.0 / (T * n_beads)
        τ = ħ * β / n_beads
    else
        τ = fixed_τ
        n_beads = Int(floor(1/(τ*T)))
    end

    #Initialsing path
    path = Path(n_beads, n_particles, n_dimensions, τ, m=m, λ = ħ^2/2/m, start_range=20.0)
    println(path.K_factor)

    if threads
        plot_on = false
    else
        plot_on = true
    end

    # Set potential function
    if pot == "Frohlich"
        potential = FrohlichPotential(α,ω,ħ) # ω is phonon frequency
    elseif pot == "Harmonic"
        potential = HarmonicPotential(ω)
    elseif pot == "MexicanHat"
        potential = MexicanHatPotential(80000.0)
    elseif pot == "Constant"
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
        observables_skip = 100
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
    println("old pos: ", path.beads[1, 1, :]) 
    PIMC(n_thermalised, n_thermalised*2, 2*n_thermalised, path, SingleMover(path), estimators, HarmonicPotential(1.0), regime, observables, adjust=true)
    if verbose
        println("Thermalisation complete")
    end
    
    println("new pos: ", path.beads[1, 1, :])
    # Comparison energy
    comparison_energy = 0.0;
    if typeof(potential) == HarmonicPotential
        comparison_energy = analyticEnergyHarmonic(potential.ω,β,ħ,n_dimensions)
    elseif typeof(potential) == FrohlichPotential

        #comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
        #comparison_energy = polaron([α], [T], [0.0]; ω=1.0, verbose = true)
        
        if length(ω) == 1
            comparison_polaron = feynmanvw(2.0, 1.0, α, ω, Inf)
        else
            comparison_polaron = feynmanvw(2.0, 1.0, [α], ω, Inf)
        end
        comparison_energy += comparison_polaron[3]
        
        #end
    end

    if !threads
        
        data = PIMC(n_steps, equilibrium_skip, observables_skip, path, mover, estimators, potential, regime, observables, adjust=true)
        
        # Storing all the different outputz
        energies = data["Energy:$(estimator)"]
        positions = data["Position:p$(particleIndex)d$(dimensionIndex)"] # Select a particular particle and dimension
        correlations = data["Correlation:$(estimator)"]
        acceptance_rates = data["Acceptance Rate:p$(particleIndex)"]
        adjuster_values = data["Adjuster Value:p$(particleIndex)"]

        # Post analysis
        variances = jackknife(energies)
        jacknife_errors = sqrt(variances[2])
        mean_energy = mean(energies)

        #=
        if n_particles == 2
            positions1 = positions
            positions2 = collect(Iterators.flatten(data["Position:p2d1"]))
            posplot = histogram([positions1, positions2])
            display(posplot)
        end
        =#

        # Output measurements and statistics
        println("Number of Beads: ", n_beads)
        println("Number of steps: ", n_steps)
        println("τ is: ", τ)
        println("Temperature: ", T)
        println("α: ", α)
        println("Mean Energy: ", mean_energy)
        println("Comparison Energy: ", comparison_energy)
        println("jackknife errors: ", jacknife_errors)
        
        # return energy, variances, mean_acceptance_rate, comparison_energy
        return data, mean_energy, jacknife_errors, comparison_energy, equilibrium_skip, observables_skip, n_beads, path
    
    else
        println("multithreading")
        aggregated_data = Dict() # Create an empty dictionary for storing data
        Mean_energy_arr = [0.0 for j in 1:thread_number]
        Error_arr = [0.0 for j in 1:thread_number]

        @threads for i in 1:thread_number
            println("i = $i on thread $(Threads.threadid())")

            data = PIMC(n_steps, equilibrium_skip, observables_skip, deepcopy(path), deepcopy(mover), estimators, potential, regime, observables, adjust=true)
        
            # Storing all the different outputz
            energies = data["Energy:$(estimator)"]
            positions = data["Position:p$(particleIndex)d$(dimensionIndex)"] # Select a particular particle and dimension
            correlations = data["Correlation:$(estimator)"]
            acceptance_rates = data["Acceptance Rate:p$(particleIndex)"]
            adjuster_values = data["Adjuster Value:p$(particleIndex)"]
            data["Comparison Energy"] = comparison_energy

            # Post analysis
            variances = jackknife(energies)
            jacknife_errors = sqrt(variances[2])
            mean_energy = mean(energies)

            # Output measurements and statistics
            aggregated_data[i] = data;
            Mean_energy_arr[i] = mean_energy;
            Error_arr[i] = jacknife_errors;
        end

        println("Number of Paths:", thread_number)
        println("Number of Beads: ", n_beads)
        println("Number of steps: ", n_steps)
        println("τ is: ", τ)
        println("Temperature: ", T)
        println("α: ", α)
        println("Mean Energy: ", mean(Mean_energy_arr))
        println("Comparison Energy: ", comparison_energy)
        println("jackknife errors: ", mean(Error_arr))

        return aggregated_data, Mean_energy_arr, Error_arr, n_beads;
    end

end
#=
function MultiModePIMC(T::Float64, ω_array::Vector{Float64}, α::Float64, n_dimensions::Float64, n_steps::Int64;
    m::Float64=1.0, n_particles::Int64=1, regime::String="Primitive", fixed_beads::Bool=false, 
    fixed_τ::Float64=0.1, n_beads::Int64=50, n_thermalised::Int64=10000, mover::String="Single", 
    potential::String="Frohlich", estimator::String="Virial", quick_steps::Bool=false, threads::Bool = false, 
    start_range::Float64 = 1.0, particleIndex::Int64 = 1, dimensionIndex::Int64 = 1,
    observable_skip_factor::Float64=0.005, equilibrium_skip_factor::Float64=0.5, version::Int64 = 1, 
    verbose::Bool = true, thread_number::Int64 = 16)
    """
    Allow array of phonon frequecies values (multi-mode)

    Initialise System Variables:
    T::Array of values
    m::Float
    n_particles::Int
    n_dimensions::Int
    start_range::Float

    Initialise potential variables:
    ω_array::Vector{Float64}
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
    β = 1 / T

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

    return 


end
=#

function RangeAlphaPIMC(T::Float64, ω::Float64, α_array::Vector{Float64}, n_dimensions::Int64, n_steps::Int64; 
    n_particles::Int64 = 1, m::Float64=1.0, regime::String="Primitive", fixed_beads::Bool=false, fixed_τ::Float64=0.1, n_beads::Int64=50, 
    n_thermalised::Int64=10000, mover::String="Single", pot::String="Frohlich", estimator::String="Virial", 
    quick_steps=false, threads::Bool = false, start_range = 1.0, particleIndex::Int64 = 1, dimensionIndex::Int64 = 1,
    observable_skip_factor::Float64=0.005, equilibrium_skip_factor::Float64=0.5, version = 1, verbose::Bool = true, 
    thread_number::Int64 = 16)
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

    @threads for i in 1:length(α_array)

        println("Simulating multiple alphas")
        α = α_array[i]

        # Set potential function
        if pot == "Frohlich"
            potential = FrohlichPotential(α,ω,ħ) # ω is phonon frequency
        elseif pot == "Harmonic"
            potential = HarmonicPotential(ω)
        elseif pot == "MexicanHat"
            potential = MexicanHatPotential(80000.0)
        elseif pot == "Constant"
            potential = ConstantPotential(10.0)
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
        PIMC(n_thermalised, n_thermalised, n_thermalised, deepcopy(path), SingleMover(deepcopy(path)), estimators, HarmonicPotential(1.0), regime, observables, adjust=true)
        if verbose
            println("Thermalisation complete")
        end

        #data = PIMC(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true, visual=false)
        data = PIMC(n_steps, equilibrium_skip, observables_skip, deepcopy(path), deepcopy(mover), estimators, potential, regime, observables, adjust=true)
        
        # Storing all the different outputz
        energies = data["Energy:$(estimator)"]
        positions = data["Position:p$(particleIndex)d$(dimensionIndex)"] # Select a particular particle and dimension

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

    return Mean_energy_arr, Error_arr, aggregated_data
end

function RangeTempPIMC(T_array::Vector{Float64}, ω::Float64, α::Float64, n_dimensions::Int64, n_steps::Int64; m::Float64=1.0, n_particles::Int64=1, 
    regime::String="Primitive", fixed_beads::Bool = false, fixed_τ::Float64 = 0.1, n_beads::Int64=50, n_thermalised::Int64=10000, 
    mover::String="Single", pot::String="Frohlich", estimator::String="Virial", quick_steps::Bool=false, threads::Bool = false, 
    start_range::Float64 = 1.0, particleIndex::Int64 = 1, dimensionIndex::Int64 = 1, observable_skip_factor=0.005, 
    equilibrium_skip_factor=0.5, version = 1, verbose::Bool = true, thread_number::Int64 = 16)
    
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
    if pot == "Frohlich"
        potential = FrohlichPotential(α,ω,ħ) # ω is phonon frequency
    elseif pot == "Harmonic"
        potential = HarmonicPotential(ω)
    elseif pot == "MexicanHat"
        potential = MexicanHatPotential(80000.0)
    elseif pot == "Constant"
        potential = ConstantPotential(10.0)
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
            mover_choice = SingleMover(path)
        elseif mover == "Displace"
            mover_choice = DisplaceMover(path)
        elseif mover == "Bisect"
            mover_choice = BisectMover(path)
        else
            println("Invalid mover")
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

        
        data = PIMC(n_steps, equilibrium_skip, observables_skip, path, mover_choice, estimators, potential, regime, observables, adjust=true)
        
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

    return Mean_energy_arr, Error_arr, aggregated_data
end

function general_Holstein_PIMC(T::Float64, ω::Float64, α::Float64, n_dimensions::Int64, n_steps::Int64; 
    m::Float64=1.0, J::Float64=1.0, n_particles::Int64 = 1, fixed_beads::Bool=false, fixed_τ::Float64 = 0.05, 
    n_beads::Int64 = 100, potential = "Holstein", estimator="Thermodynamic", 
    quick_steps=false, threads::Bool = false, particleIndex = 1, dimensionIndex = 1,
    observable_skip_factor=0.005, equilibrium_skip_factor=0.5, version = 1, 
    verbose::Bool = true, thread_number = 16)
    
    """
    Initialise System Variables:
    T
    m
    n_particles
    n_dimensions
    start_range

    Initialise potential variables:
    ω
    α

    Possible choice of potential: {"Harmonic", "Frohlich", "MexicanHat", "Constant"}
    Possible choice of estimator: {"Simple", "Virial", "Thermodynamic"}
    Possible choice of regime: {"Simple", "Primitive"}
    Possible choice of mover: {"Single", "Displace", "Bisect"}
    """
    
    #version = Int(floor(rand()*10000))

    # Path variables
    β = 1 / T
    ħ = 1.0

    # For fixed number of beads or not fixing beads
    if fixed_beads
        n_beads = n_beads
        #τ = 1.0 / (T * n_beads)
        τ = ħ * β / n_beads
    else
        τ = fixed_τ
        n_beads = Int(floor(1/(τ*T)))
    end

    # Initate path
    path = DiscretePath(n_beads, n_particles, n_dimensions, τ, m)

    if threads
        plot_on = false
    else
        plot_on = true
    end

    # Set potential function
    if potential == "Holstein"
        potential = HolsteinPotential(α, ω, ħ, J)
    end

    # Set Estimator
    if estimator == "Thermodynamic"
        estimators = [ThermodynamicEstimator()]
    else
        println("Invalid Estimator: ", estimator)
    end

    # Observables type
    observables = ["Energy", "Position"]

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

    if verbose
        println("------started--------")
        println("Temperature is: ", T)
        println("α is: ", α)
        println("n_step is: ", n_steps)
    end

    if !threads
        
        data = Holstein_PIMC(n_steps, equilibrium_skip, observables_skip, path, estimators, potential, observables)
        
        # Storing all the different outputz
        energies = data["Energy:$(estimator)"]
        positions = data["Position:p$(particleIndex)d$(dimensionIndex)"] # Select a particular particle and dimension
        #correlations = data["Correlation:$(estimator)"]

        # Post analysis
        variances = jackknife(energies)
        jacknife_errors = sqrt(variances[2])
        mean_energy = mean(energies)

        # Output measurements and statistics
        if verbose
            println("Number of Beads: ", n_beads)
            println("Number of steps: ", n_steps)
            println("τ is: ", τ)
            println("Temperature: ", T)
            println("α: ", α)
            println("Mean Energy: ", mean_energy)
            println("jackknife errors: ", jacknife_errors)
        end
        
        energy_plot = plot(energies, ylabel="Energy", xlab = "Sweeps / $observables_skip\$ n\$")
        display(energy_plot)
        # return energy, variances, mean_acceptance_rate, comparison_energy
        SaveJLDData(T, potential, version, n_beads, n_steps, data)

        return mean_energy, jacknife_errors, data
    
    else
        println("multithreading")
        aggregated_data = Dict() # Create an empty dictionary for storing data
        Mean_energy_arr = [0.0 for j in 1:thread_number]
        Error_arr = [0.0 for j in 1:thread_number]

        @threads for i in 1:thread_number
            println("i = $i on thread $(Threads.threadid())")

            data = PIMC(n_steps, equilibrium_skip, observables_skip, deepcopy(path), deepcopy(mover), estimators, potential, regime, observables, adjust=true)
        
            # Storing all the different outputz
            energies = data["Energy:$(estimator)"]
            positions = data["Position:p$(particleIndex)d$(dimensionIndex)"] # Select a particular particle and dimension
            correlations = data["Correlation:$(estimator)"]
            acceptance_rates = data["Acceptance Rate:p$(particleIndex)"]
            data["Comparison Energy"] = comparison_energy

            # Post analysis
            variances = jackknife(energies)
            jacknife_errors = sqrt(variances[2])
            mean_energy = mean(energies)

            # Output measurements and statistics
            aggregated_data[i] = data;
            Mean_energy_arr[i] = mean_energy;
            Error_arr[i] = jacknife_errors;
        end

        println("Number of Paths:", thread_number)
        println("Number of Beads: ", n_beads)
        println("Number of steps: ", n_steps)
        println("τ is: ", τ)
        println("Temperature: ", T)
        println("α: ", α)
        println("Mean Energy: ", mean(Mean_energy_arr))
        println("Comparison Energy: ", comparison_energy)
        println("jackknife errors: ", mean(Error_arr))

        return aggregated_data, Mean_energy_arr, Error_arr, n_beads;
    end

end

function SaveJLDData(T, potential, version, n_beads, n_steps, data)
    save("data_arr/$(string(typeof(potential)))/$(string(Symbol(potential)))_T$(T)_nsteps$(n_steps)_v$(version)_beads$(n_beads).jld", "data", data)
end

function SaveJLDData(T, potential, version, n_beads, n_steps, data)
    save("data_arr/$(string(typeof(potential)))/$(string(Symbol(potential)))_T$(T)_nsteps$(n_steps)_v$(version)_beads$(n_beads).jld", "data", data)
end