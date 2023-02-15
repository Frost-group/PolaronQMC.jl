# PIMC.jl
using Base.Threads


# PIMC function with multi-threading support
function PIMC(n_steps::Int, equilibrium_skip, observable_skip, path::Path, movers::Dict, observables, estimators::Array, potential::Potential, regime::Regime; adjust::Bool = true, visual::Bool = false, threads::Bool = false)
	
	# Conversion of estimators to subscripable string
	estimators_string = [string(Symbol(estimator)) for estimator in estimators]

	# Initialise independent adjusters for each mover
	adjusters = Dict()
	for mover in keys(movers)
		adjusters[mover] = path.adjusters[mover]
	end

	# Output arrays for different estimators of observables
	output_observables = Dict()
	for observable in observables
		output_observables[string(observable)] = Dict() 
		for estimator in estimators_string
			output_observables[string(observable)][estimator] = []
		end
	end

	# Output dictionary that stores width and acceptance rate
	adjuster_stats = Dict()
	for mover in keys(movers)
		adjuster_stats[mover]= Dict("Acceptance Rate" => [], "Shift Width" => [])
	end

	# Position for visuals
	visual_positions = []

	# Processes that run per step
	if threads
		@threads for step in 1:n_steps
			
			# Updating n_accepted, moving beads, and changing shift width if necessary
			for particle in rand(1:path.n_particles, path.n_particles)
				for mover_index in 1:length(movers[1])
					adjuster = path.adjusters[string(Symbol(movers[1][mover_index]))]
					if movers[2][mover_index] > rand()
						system_stats["attempted_array"][movers_string[mover_index]] += 1
						system_stats["acceptance_array"][movers_string[mover_index]] += movers[1][mover_index](path, particle, potential, regime, adjuster)
					end
				end

				# Changing shift width automatically
				if adjust 
					for adjuster in values(path.adjusters)
						update_shift_width!(adjuster)
					end
				end
			end

			# Generates observable for each cycle of "observable_skip"
			if mod(step, observable_skip) == 0 && step > equilibrium_skip

				for observable in observables
					for estimator_index in 1:length(estimators)
						append!(output_observables[string(observable)][estimators_string[estimator_index]], observable(path, potential, estimators[estimator_index]))

					end
				end

				if visual
					push!(visual_positions,copy(path.beads))
				end
			end
		end

	else
		for step in 1:n_steps
			println("step is ", step)

			# Updating n_accepted, moving beads, and changing shift width if necessary
			for particle in rand(1:path.n_particles, path.n_particles)
				#for sweep in 1:path.n_beads
				for sweep in 1:100
					for mover in keys(movers)
						if movers[mover][1] > rand()
							getfield(Main, Symbol(mover))(
								path, particle, potential, regime, adjusters[mover]
							)
						end
					end
				end
			end

			# Changing shift width automatically and save results
			if adjust
				for mover in keys(movers)
					update_shift_width!(adjusters[mover], potential)
					push!(adjuster_stats[mover]["Acceptance Rate"], adjusters[mover].acceptance_rate)
					push!(adjuster_stats[mover]["Shift Width"], adjusters[mover].value)
				end
			end

			# Generates observable for each cycle of "observable_skip"
			if mod(step, observable_skip) == 0 && step > equilibrium_skip

				for observable in observables
					for estimator_index in 1:length(estimators)
						append!(output_observables[string(observable)][estimators_string[estimator_index]], observable(path, potential, estimators[estimator_index]))

					end
				end

				# Save position in memory for animation
				if visual
					push!(visual_positions,copy(path.beads))
				end
			end
		end
	end

	return [adjuster_stats, output_observables, visual_positions]
end