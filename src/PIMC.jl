# PIMC.jl



function PIMC(n_steps::Int, equilibrium_skip, observable_skip, path::Path, movers, observables, estimators, potential::Potential, regime::Regime, adjusters; adjust = true, visual = false)
	#setting up storage of output

	#acceptance and attempt arrays for movers
		acceptance_array = Dict()
		attempted_array = Dict()
		for mover in movers[1]
			acceptance_array[string(Symbol(mover))] = 0
			attempted_array[string(Symbol(mover))] = 0
		end

	#output arrays for different estimators of observables
		output_observables = Dict()
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
				adjuster = adjusters[mover_index]
				if movers[2][mover_index] > rand()
					attempted_array[string(Symbol(movers[1][mover_index]))] += 1
					acceptance_array[string(Symbol(movers[1][mover_index]))] += movers[1][mover_index](path, particle, potential, regime, adjuster)
					
					if adjust #changing shift width automatically
						for adjuster in adjusters
							update_shift_width!(adjuster)
						end
					end
					

				end
			end
		end

		#generates observable for each cycle of "observable_skip"
		if mod(step, observable_skip) == 0 && step > equilibrium_skip
			#println("ad 2 sw ",adjusters[2].shift_width)

			for observable in observables
				for estimator in estimators
				
				observable_value = observable(path, potential, estimator) #getting value of observable
				append!(output_observables[string(observable)][string(Symbol(estimator))],observable_value)
				end
			end

		end
	end

	acceptance_ratio = Dict()

	for mover in movers[1]
		acceptance_ratio[string(Symbol(mover))] = acceptance_array[string(Symbol(mover))] / attempted_array[string(Symbol(mover))] 
	end
	
	
	return [acceptance_ratio, output_observables]

end

