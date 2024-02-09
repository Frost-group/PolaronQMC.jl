# quickruns.jl
function quickrun_frohlich(T::Float64, alpha_range, fixed_τ::Float64, thermalisation_steps::Int, steps_base::Int; thermalised::Bool = true, threads::Bool = true, verbose::Bool = true)
    #default variables for simulation
        n_particles = 1
        n_dimensions = 3
        ω = 1.0
        ħ = 1.0
        
        adjusted_beads = Int(floor(1/(fixed_τ*T)))
        

    

    #default thermalised start
    if thermalised
        thermalised_path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)
        st_potential = FrohlichPotential(alpha_range[1],ω,ħ)
        thermalised_start!(thermalised_path, st_potential, n_steps = thermalisation_steps, threads = threads, verbose = verbose)
    end

    #default simulation running
        estimators = [Virial_Estimator()]
        movers = [[Bisect!],[1.0]]
        regime = Primitive_Regime()
        observables = [Energy]

        #output

        observables_range_L = Dict()
        for estimator in estimators
            observables_range_L[string(Symbol(estimator))] = []
        end

        errors_range_L = Dict()
        for estimator in estimators
            errors_range_L[string(Symbol(estimator))] = []
        end

        #loop for running of simulation
        for L in alpha_range
        
            α = L
            n_steps = Int(steps_base*L)
            potential = FrohlichPotential(α,ω,ħ)
            equilibrium_skip = 0.7*n_steps
            observables_skip = 0.03*n_steps
            
            if thermalised
                path = deepcopy(thermalised_path)
            else
                path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)
            end

            pimc = PIMC(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true, threads=threads)
    
            output_observables = pimc[2]
            totalEnergy = output_observables["Energy"]
    
            for estimator in estimators
                estimator_energy = energy[string(Symbol(estimator))]
                
                #energy
                append!(observables_range_L[string(Symbol(estimator))], mean(estimator_energy))
    
                #errors
                error = sqrt(jackknife(estimator_energy)[2])
                append!(errors_range_L[string(Symbol(estimator))], error)
    
            end
        if verbose
            println("α = $L PIMC Complete.")
        end
        end

    return [observables_range_L[string(Symbol(estimators[1]))], errors_range_L[string(Symbol(estimators[1]))]]

end

