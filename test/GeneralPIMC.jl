using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility
using LaTeXStrings
using DelimitedFiles
using JLD
using Base.Threads

function generalPIMC(T, m, ω, α, n_particles, n_dimensions, regime, fixed_beads, fixed_τ, n_beads, n_steps, 
    n_thermalised, mover, potential, estimator, quick_steps=false, threads::Bool = false, start_range = 1.0, particleIndex = 1, dimensionIndex = 1,
    observable_skip_factor=0.005, equilibrium_skip_factor=0.5, version = 1, verbose::Bool = true, thread_number = 16)
    
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

    # Comparison energy
    if typeof(potential) == HarmonicPotential
        comparison_energy = analyticEnergyHarmonic(potential.ω,β,ħ,n_dimensions)
    elseif typeof(potential) == FrohlichPotential
        comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
        comparison_energy = comparison_polaron.F
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
        return mean_energy, jacknife_errors, comparison_energy, energies, positions, correlations, acceptance_rates, adjuster_values, equilibrium_skip, observables_skip, version, potential, n_beads, data
    
    else
        println("multithreading")
        aggregated_data = Dict() # Create an empty dictionary for storing data
        Mean_energy_arr = [0.0 for j in 1:thread_number]
        Error_arr = [0.0 for j in 1:thread_number]

        @threads for i in 1:thread_number
            println("i = $i on thread $(Threads.threadid())")

            data = PIMC(n_steps, equilibrium_skip, observables_skip, path, mover, estimators, potential, regime, observables, adjust=true)
        
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

        return aggregated_data, Mean_energy_arr, Error_arr;
    end

end

function SaveJLDData(T, potential, version, n_beads, n_steps, data)
    save("data_arr/$(string(typeof(potential)))/$(string(Symbol(potential)))_T$(T)_nsteps$(n_steps)_v$(version)_beads$(n_beads).jld", "data", data)
end
#=

=#