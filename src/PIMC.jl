# PIMC.jl



function PIMC(n_steps::Int, equilibrium_skip, observable_skip, path::Path, movers, observables, estimators, potential::Potential, regime::Regime; adjust = true, visual = false)
	#setting up storage of output

	# information about the running of the system
		system_stats = Dict()
		system_stats["acceptance_array"] = Dict()
		system_stats["attempted_array"] = Dict()

	#acceptance and attempt arrays for movers

		for mover in movers[1]
			system_stats["acceptance_array"][string(Symbol(mover))] = 0
			system_stats["attempted_array"][string(Symbol(mover))] = 0
		end

	#output arrays for different estimators of observables
		output_observables = Dict()
		#observables_counter = 0 #used in weighting of the rolling averages of observables
		#generating lists for output
		for observable in observables
			output_observables[string(observable)] = Dict() 
			for estimator in estimators
				output_observables[string(observable)][string(Symbol(estimator))] = []
			end
		end

		

	
	for step in 1:n_steps
		if step == 0.5*n_steps
			println("50% complete")
		end

		#updating n_accepted, moving beads, and changing shift width if necessary
		for particle in rand(1:path.n_particles, path.n_particles)
			for mover_index in 1:length(movers[1])
				adjuster = path.adjusters[string(Symbol(movers[1][mover_index]))]
				if movers[2][mover_index] > rand()
					system_stats["attempted_array"][string(Symbol(movers[1][mover_index]))] += 1
					system_stats["acceptance_array"][string(Symbol(movers[1][mover_index]))] += movers[1][mover_index](path, particle, potential, regime, adjuster)
					

				end
			end

			if adjust #changing shift width automatically
				for adjuster in values(path.adjusters)
					update_shift_width!(adjuster)
				end
			end
		end

		#generates observable for each cycle of "observable_skip"
		if mod(step, observable_skip) == 0 && step > equilibrium_skip
			#observables_counter += 1

			for observable in observables
				for estimator in estimators
					append!(output_observables[string(observable)][string(Symbol(estimator))], observable(path, potential, estimator))

				
				#output_observables[string(observable)][string(Symbol(estimator))] = (((observables_counter-1) * output_observables[string(observable)][string(Symbol(estimator))]) + observable(path, potential, estimator))/observables_counter  #rolling weighted mean of observable
				
				end
			end

		end
	end

	system_stats["acceptance_ratio"] = Dict()

	for mover in movers[1]
		system_stats["acceptance_ratio"][string(Symbol(mover))] = system_stats["acceptance_array"][string(Symbol(mover))] / system_stats["attempted_array"][string(Symbol(mover))] 
	end
	
	
	return [system_stats["acceptance_ratio"], output_observables]

end

function primitive_action(path::Path, bead::Int, particle::Int, potentials::Array{Potential})
    primitive_action = kinetic_action(path, bead, particle)
    for potential in potentials
        primitive_action += potential_action(path, bead, particle, potential)
    end
    return primitive_action
end

function PIMC(n_steps::Int, path::Path, movers, observables, potentials::Union{Potential, Array{Potential}})
	
	observable_skip = 0.001 * n_steps
	equilibrium_skip = 0.2 * n_steps
	# equilibrium_skip = 0
	
	n_accepted = Dict(string(Symbol(mover)) => 0 for mover in movers)
	observable_traces = Dict(string(Symbol(observable)) => [] for observable in observables)
	
	path_trace = []
	for step in 1:n_steps

		for mover in movers
			for particle in rand(1:path.n_particles, path.n_particles)
				n_accepted[string(Symbol(mover))] += mover(path, particle, potentials)
			end
		end

		if mod(step, observable_skip) == 0 && step > equilibrium_skip
			for observable in observables
				append!(observable_traces[string(Symbol(observable))], [observable(path, potentials)])
			end
			append!(path_trace, [path.beads])
		end
		
	end

	acceptance_ratio = Dict(string(Symbol(mover)) => 1.0 * n_accepted[string(Symbol(mover))] / (n_steps * path.n_particles) for mover in movers)

	return acceptance_ratio, path_trace, observable_traces
end

