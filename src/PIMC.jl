# PIMC.jl
using Base.Threads


# PIMC function with multi-threading support
function PIMC(n_steps::Int, equilibrium_skip, observable_skip, path::Path, movers::Array, observables, estimators::Array, potential::Potential, regime::Regime; adjust::Bool = true, visual::Bool = false, threads::Bool = false)
	
	# Setting up storage of output

	# Conversion of objects to strings to dictionaries
	movers_string = [string(Symbol(mover)) for mover in movers[1]]
	estimators_string = []
	for estimator in estimators
		push!(estimators_string, string(Symbol(estimator)))
	end

	# Information about the running of the system
	system_stats = Dict()
	system_stats["acceptance_array"] = Dict()
	system_stats["attempted_array"] = Dict()

	# Acceptance and attempt arrays for movers
	for mover in movers_string
		system_stats["acceptance_array"][mover] = 0
		system_stats["attempted_array"][mover] = 0
	end

	# Output arrays for different estimators of observables

	output_observables = Dict()
	# Generating lists for output
	for observable in observables
		output_observables[string(observable)] = Dict() 
		for estimator in estimators_string
			output_observables[string(observable)][estimator] = []
		end
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
			
			# Updating n_accepted, moving beads, and changing shift width if necessary
			for particle in rand(1:path.n_particles, path.n_particles)
				for mover_index in 1:length(movers[1])
					adjuster = path.adjusters[movers_string[mover_index]]
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
	end

	system_stats["acceptance_ratio"] = Dict()

	for mover in movers_string
		system_stats["acceptance_ratio"][mover] = system_stats["acceptance_array"][mover] / system_stats["attempted_array"][mover] 
	end

	return [system_stats["acceptance_ratio"], output_observables, visual_positions]
end