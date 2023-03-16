# PIMC.jl
using Base.Threads
using StaticArrays
using PyCall


# PIMC function with multi-threading support
function PIMC(n_steps::Int, equilibrium_skip, observable_skip, path::Path, mover::Mover, estimators::Array, potential::Potential, regime::Regime, observables::Array; adjust::Bool = true, threads::Bool = false)

	# Conversion of estimators to subscripable string
	estimators_string = [split(string(Symbol(estimator)), "Estimator()")[1] for estimator in estimators]
	observables_set = Set(observables)

	# Dictionary to store all PIMC data
	data = Dict()

	# If we consider Frohlich as bounded, we need to consider the initial bead and final bead differently
	if typeof(regime) == BoundRegime
		n_sweep = path.n_beads+1
	else
		n_sweep = path.n_beads
	end

	# If we are using Bisect function, then we can reduce the sweep. The power of 2 can be changed, but also need to change the moves.jl
	if typeof(mover) == BisectMover
		n_sweep = Int(floor(path.n_beads/2^5))
	end

	# Create data structures for energies
	if "Energy" in observables_set
		for estimator in estimators_string
			data["Energy:$(estimator)"] = []
		end
	end

	# Create data structures for positions
	if "Position" in observables_set
		for particle in 1:path.n_particles
			for dimension in 1:path.n_dimensions
				data["Position:p$(particle)d$(dimension)"] = []
			end
		end
	end

	if "Correlation" in observables_set
		for estimator in estimators_string
			data["Correlation:$(estimator)"] = []
		end
	end

	# Create data structures for adjuster stats
	if adjust
		for particle in 1:path.n_particles
			data["Acceptance Rate:p$(particle)"] = fill(NaN, n_steps)
			data["Adjuster Value:p$(particle)"] = fill(NaN, n_steps)
		end
	end

	# Processes that run per step
	if threads
		@threads for step in 1:n_steps

			# Updating n_accepted, moving beads, and changing shift width if necessary
			for particle in rand(1:path.n_particles, path.n_particles)
				for sweep in 1:path.n_beads
					moveBead(mover, path, particle, potential, regime)
				end
			
				# Changing shift width automatically and save results
				if adjust
					updateAdjuster(mover.adjusters[particle], potential)
					data["Acceptance Rate"][particle, step] = mover.adjusters[particle].acceptance_rate
					data["Adjuster Value"][particle, step] = mover.adjusters[particle].value
				end
			end

			# Generates observable for each cycle of "observable_skip"
			if mod(step, observable_skip) == 0 && step > equilibrium_skip

				# Add energy of system for step to data
				if energies
					for (index, estimator) in enumerate(estimators_string)
						data["Energy"][estimator][step] = Energy(path, potential, estimators[index])
					end
				end

				# Add positions for step to data for each particle
				if positions
					for particle in 1:path.n_particles
						push!(data["Position"][particle], copy(path.beads[:, particle, :]))
					end
				end
			end
		end
	else
		for step in 1:n_steps
			println("step is: ", step)

			# Updating n_accepted, moving beads, and changing shift width if necessary
			for particle in rand(1:path.n_particles, path.n_particles)
				for sweep in 1:n_sweep
					moveBead(mover, path, particle, potential, regime)
				end
			
				# Changing shift width automatically and save results
				if adjust
					updateAdjuster(mover.adjusters[particle], potential)
					data["Acceptance Rate:p$(particle)"][step] = mover.adjusters[particle].acceptance_rate
					data["Adjuster Value:p$(particle)"][step] = mover.adjusters[particle].value
				end
			end

			# Generates observable for each cycle of "observable_skip"
			if mod(step, observable_skip) == 0 && step > equilibrium_skip

				# Add energy of system for step to data
				if "Energy" in observables_set
					for (index, estimator) in enumerate(estimators_string)
						push!(data["Energy:$(estimator)"], energy(path, potential, estimators[index]))
					end
				end

				# Add positions for step to data for each particle
				if "Position" in observables_set
					for particle in 1:path.n_particles
						for dimension in 1:path.n_dimensions
							push!(data["Position:p$(particle)d$(dimension)"], copy(path.beads[:, particle, dimension]))
						end
					end
				end

				# Add correlation to data fro each estimator
				if "Correlation" in observables_set
					for (estimator, estimator_string) in zip(estimators, estimators_string)
						push!(data["Correlation:$(estimator_string)"], correlation(path, potential, estimator))
					end
				end
			end
		end
	end

	return data
end