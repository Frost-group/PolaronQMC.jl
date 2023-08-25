# PIMC.jl
using Base.Threads
using StaticArrays
using PyCall

using TimerOutputs
const tmr = TimerOutput();

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
	#tmr = TimerOutput();
	# Conversion of estimators to subscripable string
	estimators_string = [split(string(Symbol(estimator)), "Estimator()")[1] for estimator in estimators]
	observables_set = Set(observables)

	# Dictionary to store all PIMC data; The data are store in format of data["String"] 
	data = Dict()

	n_sweep = path.n_beads

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

	# Preallocation of memories
	action_arr = [0.0, 0.0] #index 1: Old action; index 2: New Action
	shift = zeros(path.n_dimensions) # Shifting sizes
	store_diff = zeros(path.n_dimensions) # used to store the "temporary array" when subtracting two arrays
	
	# Preallocation cosh entries as it remains constant
	coshfn = zeros(path.n_beads, path.n_beads, length(potential.ω))
	if typeof(potential) == FrohlichPotential
		for k in 1:length(potential.ω)
			fac = potential.ħ * path.τ * path.n_beads * potential.ω[k]; 
			for i in 1:path.n_beads, j in i+1:path.n_beads
				coshfn[j, i, k] = cosh(fac * ((j-i)/path.n_beads-0.5))
			end
		end
	end


	# Processes that run per step (Only single threaded kept as multi-threading at this layer create problem)
	for step in 1:n_steps

		# Feedback for the progress of the PIMC simulation
		if mod(step, 2000) == 0
			println("step is: ", step)
		end

		if mod(step, 100) == 0
			if typeof(potential) == FrohlichPotential
				# Want to shift it back such that it centres at origin
				recentralise(path, verbose=false)
			end
		end

		# Updating n_accepted, moving beads, and changing shift width if necessary
		for particle in rand(1:path.n_particles, path.n_particles)
			#@timeit tmr "Moving" 
			for sweep in 1:n_sweep
				moveBead!(mover, path, particle, potential, regime, action_arr, shift, store_diff, coshfn) # Moving beads a total of n_sweep times
			end
		
			# Changing shift width automatically and save results
			if adjust
				#@timeit tmr "Adjusting" 
				updateAdjuster(mover.adjusters[particle], path)
				data["Acceptance Rate:p$(particle)"][step] = mover.adjusters[particle].acceptance_rate
				data["Adjuster Value:p$(particle)"][step] = mover.adjusters[particle].value
			end
		end

		# Generates observable for each cycle of "observable_skip"
		if mod(step, observable_skip) == 0 && step > equilibrium_skip

			# Add energy of system for step to data
			#@timeit tmr "Cal Energy" 
			if "Energy" in observables_set
				for (index, estimator) in enumerate(estimators_string)
					push!(data["Energy:$(estimator)"], energy(path, potential, estimators[index], store_diff, coshfn))
				end
			end

			# Add positions for step to data for each particle
			#@timeit tmr "Cal Position" 
			if "Position" in observables_set
				for particle in 1:path.n_particles
					for dimension in 1:path.n_dimensions
						push!(data["Position:p$(particle)d$(dimension)"], collect(Iterators.flatten(path.beads[:, particle, dimension])))
					end
				end
			end

			# Add correlation to data fro each estimator
			#@timeit tmr "Cal Correlation" 
			if "Correlation" in observables_set
				for (estimator, estimator_string) in zip(estimators, estimators_string)
					push!(data["Correlation:$(estimator_string)"], correlation(path, potential, estimator))
				end
			end
		end
	end
	
	#println(tmr)
	reset_timer!(tmr::TimerOutput)
	return data
end

function Holstein_PIMC(n_steps::Int, equilibrium_skip, observable_skip, path::DiscretePath, estimators::Array, potential::Potential, observables::Array)
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
	estimators_string = [split(string(Symbol(estimator)), "Estimator()")[1] for estimator in estimators]
	observables_set = Set(observables)

	# Dictionary to store all PIMC data; The data are store in format of data["String"] 
	data = Dict()

	n_sweep = 1

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

	F_l = zeros(path.n_beads) # should start with zero distance
	for l in 1:path.n_beads
		F_l[l] = path.τ^3 * potential.λ^2 / (4 * path.n_beads) * sum((cos(2π*(k-1)*(l-1)/path.n_beads)/(1 - cos(2π*(k-1)/path.n_beads) + path.τ^2 * potential.ω^2 / 2)) for k in 1:path.n_beads)
	end

    DF_l = zeros(path.n_beads) # should start with zero distance
	for l in 1:path.n_beads
		DF_l[l] = path.τ^2 * potential.λ^2 / (4 * path.n_beads) * 
                    (3 * sum((cos(2π*(k-1)*(l-1)/path.n_beads)/(1 - cos(2π*(k-1)/path.n_beads) + path.τ^2 * potential.ω^2 / 2)) for k in 1:path.n_beads) - 
                    path.τ^2 * potential.ω^2 * sum((cos(2π*(k-1)*(l-1)/path.n_beads)/(1 - cos(2π*(k-1)/path.n_beads) + path.τ^2 * potential.ω^2 / 2)^2) for k in 1:path.n_beads))
		
	end


    particle = 1
	# Processes that run per step (Only single threaded kept as multi-threading at this layer create problem)
	for step in 1:n_steps
		#println("step is: ", step)
		# Feedback for the progress of the PIMC simulation
		if mod(step, 200000) == 0
			#plot(path.beads[:, particle, 1], 1:path.n_beads)
			println("step is: ", step)
		end



        for sweep in 1:n_sweep
            moveBead!(path, particle, potential, F_l) # Moving beads a total of n_sweep times
        end

		# Generates observable for each cycle of "observable_skip"
		if mod(step, observable_skip) == 0 && step > equilibrium_skip

			# Add energy of system for step to data
			if "Energy" in observables_set
				for (index, estimator) in enumerate(estimators_string)
					push!(data["Energy:$(estimator)"], energy(path, potential, estimators[index], DF_l))
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