# PIMC.jl

function relabel_beads!(path::Path)
	rand_slice = rand(1:path.n_beads)
	slices = vcat(rand_slice:path.n_beads, 1:rand_slice-1)
	path.beads = path.beads[slices, :, :]
end

function PIMC(n_steps::Int, equilibrium_skip, observable_skip, path::Path, movers, observables, estimator::Estimator, potential::Potential, regime::Regime; visual = false)
	
	#observable_skip = 0.01 * n_steps
	#equilibrium_skip = 0.1 * n_steps
	
	acceptance_array = Dict()
	attempted_array = Dict()
	for mover in movers[1]
		acceptance_array[string(Symbol(mover))] = 0
		attempted_array[string(Symbol(mover))] = 0
	end


	output_observables = Dict()

	#generating lists for output
	for observable in observables
		output_observables[string(observable)] = [] 
	end
	

	path_trace = []
	for step in 1:n_steps
		if step == 0.5*n_steps
			println("50% complete")
		end

		#updating n_accepted with yes or no, and moving beads
		for particle in rand(1:path.n_particles, path.n_particles)
			for mover_index in 1:length(movers[1])
				if movers[2][mover_index] > rand()
					attempted_array[string(Symbol(movers[1][mover_index]))] += 1
					acceptance_array[string(Symbol(movers[1][mover_index]))] += movers[1][mover_index](path, particle, potential, regime)
				end
			end
		end

		#generates observable for each cycle of "observable_skip"
		if mod(step, observable_skip) == 0 && step > equilibrium_skip

			for observable in observables
				
				observable_value = observable(path, potential, estimator) #getting value of observable
				append!(output_observables[string(observable)],observable_value)

			end

		end
	end

	acceptance_ratio = Dict()

	for mover in movers[1]
		acceptance_ratio[string(Symbol(mover))] = acceptance_array[string(Symbol(mover))] / attempted_array[string(Symbol(mover))] 
	end
	
	
	return [attempted_array, acceptance_ratio, output_observables]

end

