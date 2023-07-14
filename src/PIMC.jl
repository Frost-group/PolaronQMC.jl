# PIMC.jl
using Base.Threads
using StaticArrays
using PyCall


# PIMC function with multi-threading support
function PIMC(n_steps::Int, equilibrium_skip, observable_skip, path::Path, mover::Mover, estimators::Array, potential::Potential, regime::Regime, observables::Array; adjust::Bool = true, threads::Bool = false, max_level::Int = 1)
	"""
	Essential parameter
		n_steps: Total number of PIMC steps that we want to used
		equilibrium_skip: Number of steps that we do not want to collect any measurements
		observable_skip: Number of steps while performing measurements
		path: The path containing all particles' and beads' location, and some other parameter, e.g. mass, λ
		mover: Type of mover within the set {Single, Bisect, Displace}
		estimators: Energy estimator within the set {Simple, Thermodynamic, Virial}
		potential: Type of potential
		regime: Type of regime for the PIMC simulation {primitive, simple}
		observables: Type of measurements that will be extracted, e.g. Position, Correlation function, Energy
	
	Optional parameter
		adjust (Bool, default True): Decide whether to implement adjustor for a faster convergence 
		threads (Bool, default False): Choice of Multithreading. Turned off strictly due to algorithmic clashes which lead to incorrect simulation
		max_level (Int, default 1): Set the max bisecting level for the case of BisectMover (won't affect SingleMover and DisplaceMover)
	
	"""

	# Conversion of estimators to subscripable string
	estimators_string = [split(string(Symbol(estimator)), "Estimator()")[1] for estimator in estimators]
	observables_set = Set(observables)

	# Dictionary to store all PIMC data; The data are store in format of data["String"] 
	data = Dict()

	# Steps to set the number of sweep per step. If single, then on average per step we will move each bead once to avoid possible bias
	# If we are using Bisect function, then we can reduce the sweep. The power of 2 can be changed, but also need to change "moves.jl"
	if typeof(mover) == BisectMover
		n_sweep = Int(floor(path.n_beads/(2^max_level-1)))
	else
		n_sweep = path.n_beads
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

	# Create data structures for Positional Correlation
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

	# Processes that run per step (Only single threaded kept as multi-threading at this layer create problem)
	for step in 1:n_steps

		# Feedback for the progress of the PIMC simulation
		if mod(step, 2000) == 0
			println("step is: ", step)
		end

		# Updating n_accepted, moving beads, and changing shift width if necessary
		for particle in rand(1:path.n_particles, path.n_particles)
			for sweep in 1:n_sweep
				moveBead(mover, path, particle, potential, regime, max_level) # Moving beads a total of n_sweep times
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
						push!(data["Position:p$(particle)d$(dimension)"], collect(Iterators.flatten(path.beads[:, particle, dimension])))
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
	
	return data
end