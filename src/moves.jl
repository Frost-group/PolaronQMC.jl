# moves.jl
using Distributions

function moveBead(mover::SingleMover, path::Path, particle::Int, potential::FrohlichPotential, regime::BoundRegime)
	
	"""
	Single!(path::Path, particle::Int, potential::Potential, regime::Regime, adjuster::Adjuster)

	Move a single imaginary-time timeslice (bead) on a single particle, or subset of particles, per Monte-Carlo iteration.

	Arguments
	- `path::Path`: collection of all particle imaginary-time paths.
	- `particle::Union{Int, Vector{Int}}`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
	- `potentials::Union{Potential, Array{Potential}}`: list of potentials active in the system. Can just be a single potential.

	See also [`Path`](@ref).
	"""

    bead = rand(1:path.n_beads+1)	
	width = mover.adjusters[particle].value
	shift = rand(path.n_dimensions) * width .* rand([-1,1],path.n_dimensions)			# Linear random displacement of bead.
	
	# Attempt one Single move
	mover.adjusters[particle].attempt_counter += 1
	
	#if bead != (path.n_beads+1) || bead != 1
    # We just need to look at the beads +- 1 unit from m
    # CHECK: Non local potential? Coulombic?
	old_action = 
		kineticAction(path, bead-1, bead, particle, regime) +		# Link bead-1 to bead
		kineticAction(path, bead, bead+1, particle, regime) +		# Link bead to bead+1
		potentialAction(path, bead, particle, potential, regime)	# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.

	path.beads[mod1(bead, path.n_beads), particle, :] += shift

	new_action =
		kineticAction(path, bead-1, bead, particle, regime) +		# Link bead-1 to bead
		kineticAction(path, bead, bead+1, particle, regime) +		# Link bead to bead+1
		potentialAction(path, bead, particle, potential, regime)	# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.
	# Metropolis algorithm. 
	# Accept if bead displacement decreases the action, otherwise accept with probability exp(-ΔAction).

	if new_action - old_action <= 0.0 || rand() <= exp(-(new_action - old_action))
		mover.adjusters[particle].success_counter += 1 # Updating counter for adjustment of shift width
		return true
	else
		path.beads[mod1(bead, path.n_beads), particle, :] -= shift
		return false
	end
end


function moveBead(mover::SingleMover, path::Path, particle::Int, potential::Potential, regime::Regime)
	
	"""
	Single!(path::Path, particle::Int, potential::Potential, regime::Regime, adjuster::Adjuster)

	Move a single imaginary-time timeslice (bead) on a single particle, or subset of particles, per Monte-Carlo iteration.

	Arguments
	- `path::Path`: collection of all particle imaginary-time paths.
	- `particle::Union{Int, Vector{Int}}`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
	- `potentials::Union{Potential, Array{Potential}}`: list of potentials active in the system. Can just be a single potential.

	See also [`Path`](@ref).
	"""

    bead = rand(1:path.n_beads)	
	width = mover.adjusters[particle].value
	shift = rand(path.n_dimensions) * width .* rand([-1,1],path.n_dimensions)			# Linear random displacement of bead.

	# Attempt one Single move
	mover.adjusters[particle].attempt_counter += 1

    # We just need to look at the beads +- 1 unit from m
    # CHECK: Non local potential? Coulombic?
    old_action = 
		kineticAction(path, bead-1, bead, particle, regime) +		# Link bead-1 to bead
		kineticAction(path, bead, bead+1, particle, regime) +		# Link bead to bead+1
		potentialAction(path, bead, particle, potential, regime)	# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.

	path.beads[bead, particle, :] += shift

    new_action =
		kineticAction(path, bead-1, bead, particle, regime) +		# Link bead-1 to bead
		kineticAction(path, bead, bead+1, particle, regime) +		# Link bead to bead+1
		potentialAction(path, bead, particle, potential, regime)	# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.

	# Metropolis algorithm. 
	# Accept if bead displacement decreases the action, otherwise accept with probability exp(-ΔAction).

	if new_action - old_action <= 0.0 || rand() <= exp(-(new_action - old_action))
		mover.adjusters[particle].success_counter += 1 # Updating counter for adjustment of shift width
		return true
	else
		path.beads[mod1(bead, path.n_beads), particle, :] -= shift
		return false
	end
end


function moveBead(mover::DisplaceMover, particle::Int, potential::Potential, regime::Regime)

	"""
	Displace!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})

	Move the entire imaginary-time timeslice (all the beads) on a single particle, or subset of particles, by the same vector per Monte-Carlo iteration.

	Arguments
	- `path::Path`: collection of all particle imaginary-time paths.
	- `particle::Union{Int, Vector{Int}}`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
	- `potentials::Union{Potential, Array{Potential}}`: list of potentials active in the system. Can just be a single potential.

	See also [`Path`](@ref).
	"""

	width = mover.adjusters[particle].value
	shift = rand(path.n_dimensions) * width .* rand([-1,1],path.n_dimensions)

	# Attempt one Single move
	mover.adjusters[particle].attempt_counter += 1

	# Save configuration of paths. Return to this configuration if Metropolis rejects the new configuration.
	old_beads = copy(path.beads[:, particle, :])

    # Evaluate old action. 
	# The kinetic action is unchanged because the relative bead positions are preserved, so we only need to evaluate the potential action at each bead.
    old_action = sum(potentialAction(path, bead, particle, potential, regime) for bead in 1:path.n_beads)
	
	# Displace every bead on each specified particle by the same random vector.
	for bead in 1:path.n_beads
		path.beads[bead, particle, :] += shift
	end

    new_action = sum(potentialAction(path, bead, particle, potential, regime) for bead in 1:path.n_beads)

	if new_action - old_action <= 0.0
		mover.adjusters[particle].success_counter += 1
		return true
	elseif rand() <= exp(-(new_action - old_action))
		mover.adjusters[particle].success_counter += 1
		return true
	else
		path.beads[:, particle, :] = old_beads
		return false
	end
end

function moveBead(mover::BisectMover, path::Path, particle::Int, potential::FrohlichPotential, regime::Regime)

	# Segment length to perform Bisect
	max_level = 5
	segment_length = Int((2^max_level))

	start_bead = rand(1:path.n_beads)
	#old_beads = deepcopy(path.beads[:,particle, :])
	old_beads = path.beads[:,particle, :]

	# Call parameters once
	τλ = path.τ * path.λ
	n_beads = path.n_beads
	n_dim = path.n_dimensions
	beadrange = collect(start_bead+1:start_bead+segment_length-1)
	for i in 1:segment_length-1
		beadrange[i] = mod1(beadrange[i], n_beads)
	end

	# Attempt one Single move
	mover.adjusters[particle].attempt_counter += 1

	total_old_action = kineticAction(path, start_bead, start_bead+1, particle, regime)
	for bead in start_bead+1:start_bead+segment_length-1
		total_old_action += kineticAction(path, bead, bead+1, particle, regime) +
							bisectPotentialAction(path, bead, beadrange, particle, potential, regime)
	end


	for level in max_level:-1:1
		level_fac = 2^(level-1)

		# Number of Beads in Level
		beads_level = 2^(max_level - level)

		for k in 1:beads_level

			# Find bead of interest
			if k == 1
				#bead = Int(start_bead + (2^(level-1) * k))
				bead = Int(start_bead + (level_fac * k))
			else
				#bead = Int(start_bead + (2^(level-1) + (2^(level) * (k-1))))
				bead = Int(start_bead + (level_fac + (2^(level) * (k-1))))
			end

			# Find beads of which to find midpoint
			r0, r1 = bead - level_fac, bead + level_fac

			# Find normally distributed shift
			#width = sqrt(level_fac * τλ)
			#shift = width .* rand(Distributions.Normal(0, 1), n_dim)

			# Perform move and calculate change to action
			#midpoint = 0.5 * (path.beads[mod1(r0, n_beads), particle, :] + path.beads[mod1(r1, n_beads), particle, :])
			path.beads[mod1(bead, n_beads), particle, :] = 0.5 * (path.beads[mod1(r0, n_beads), particle, :] + path.beads[mod1(r1, n_beads), particle, :]) + sqrt(level_fac * τλ) .* rand(Distributions.Normal(0, 1), n_dim)
		end
	end

	total_new_action = kineticAction(path, start_bead, start_bead+1, particle, regime)
	for bead in start_bead+1:start_bead+segment_length-1
		total_new_action += kineticAction(path, bead, bead+1, particle, regime) +
							bisectPotentialAction(path, bead, beadrange, particle, potential, regime)
	end

	if total_new_action - total_old_action < 0.0 || rand() < exp(-(total_new_action - total_old_action))
		mover.adjusters[particle].success_counter += 1
		return true

	#elseif rand() < exp(-(total_new_action - total_old_action))
		#mover.adjusters[particle].success_counter += 1
		#return true
		
	else
		path.beads[:, particle, :] = old_beads
		return false
	end
end

function moveBead(mover::BisectMover, path::Path, particle::Int, potential::Potential, regime::Regime)

	# Segment length to perform Bisect
	max_level = 5
	segment_length = Int((2^max_level) + 1)

	start_bead = rand(1:path.n_beads)
	old_beads = deepcopy(path.beads[:,particle, :])

	# Call parameters once
	τλ = path.τ * path.λ
	n_beads = path.n_beads
	n_dim = path.n_dimensions

	# Attempt one Single move
	mover.adjusters[particle].attempt_counter += 1

	total_old_action = 0.0
	for bead in start_bead:start_bead+segment_length
		total_old_action += potentialAction(path, bead, particle, potential, regime)
	end

	for level in max_level:-1:1
		level_fac = 2^(level-1)

		# Number of Beads in Level
		beads_level = 2^(max_level - level)

		for k in 1:beads_level

			# Find bead of interest
			if k == 1
				#bead = Int(start_bead + (2^(level-1) * k))
				bead = Int(start_bead + (level_fac * k))
			else
				#bead = Int(start_bead + (2^(level-1) + (2^(level) * (k-1))))
				bead = Int(start_bead + (level_fac + (2^(level) * (k-1))))
			end

			# Find beads of which to find midpoint
			r0, r1 = bead - level_fac, bead + level_fac

			# Find normally distributed shift
			width = sqrt(level_fac * τλ) * mover.adjusters[particle].value
			shift = width .* rand(Distributions.Normal(0, 1), n_dim)

			# Perform move and calculate change to action
			midpoint = 0.5 * (path.beads[mod1(r0, n_beads), particle, :] + path.beads[mod1(r1, n_beads), particle, :])
			path.beads[mod1(bead, n_beads), particle, :] = midpoint + shift
		end
	end

	total_new_action = 0.0
	for bead in start_bead:start_bead+segment_length
		total_new_action += potentialAction(path, bead, particle, potential, regime)
	end

	if total_new_action - total_old_action < 0.0
		mover.adjusters[particle].success_counter += 1
		return true

	elseif rand() < exp(-(total_new_action - total_old_action))
		mover.adjusters[particle].success_counter += 1
		return true
		
	else
		path.beads[:, particle, :] = old_beads
		return false
	end
end