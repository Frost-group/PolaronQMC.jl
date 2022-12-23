# moves.jl
using Distributions

function Single!(path::Path, particle::Int, potential::Potential, regime::Regime , adjuster::Adjuster)
	
	"""
	Single!(path::Path, particle::Int, potential::Potential, regime::Regime, adjuster::Adjuster)

	Move a single imaginary-time timeslice (bead) on a single particle, or subset of particles, per Monte-Carlo iteration.

	Arguments
	- `path::Path`: collection of all particle imaginary-time paths.
	- `particle::Union{Int, Vector{Int}}`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
	- `potentials::Union{Potential, Array{Potential}}`: list of potentials active in the system. Can just be a single potential.

	See also [`Path`](@ref).
	"""

    bead = rand(1:path.n_beads)							# Pick a random bead.
	width = adjuster.shift_width
	#width = 2
	shift = rand(path.n_dimensions) * width .* rand([-1,1],path.n_dimensions)			# Linear random displacement of bead.

	# Attempt one Single move
	adjuster.attempt_counter += 1

    # We just need to look at the beads +- 1 unit from m
    # CHECK: Non local potential? Coulombic?
    old_action = 
		kinetic_action(path, bead-1, bead, particle, regime) +		# Link bead-1 to bead
		kinetic_action(path, bead, bead+1, particle, regime) +		# Link bead to bead+1
		potential_action(path, bead, particle, potential, regime)	# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.

	path.beads[bead, particle, :] += shift

    new_action =
		kinetic_action(path, bead-1, bead, particle, regime) +		# Link bead-1 to bead
		kinetic_action(path, bead, bead+1, particle, regime) +		# Link bead to bead+1
		potential_action(path, bead, particle, potential, regime)	# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.

	# Metropolis algorithm. 
	# Accept if bead displacement decreases the action, otherwise accept with probability exp(-ΔAction).

	if new_action - old_action <= 0.0 || rand() <= exp(-(new_action - old_action))
		adjuster.success_counter += 1 # Updating counter for adjustment of shift width
		return true
	else
		path.beads[bead, particle, :] -= shift
		return false
	end
end


function Displace!(path::Path, particle::Int, potential::Potential, regime::Regime, adjuster::Adjuster)

	"""
	Displace!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})

	Move the entire imaginary-time timeslice (all the beads) on a single particle, or subset of particles, by the same vector per Monte-Carlo iteration.

	Arguments
	- `path::Path`: collection of all particle imaginary-time paths.
	- `particle::Union{Int, Vector{Int}}`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
	- `potentials::Union{Potential, Array{Potential}}`: list of potentials active in the system. Can just be a single potential.

	See also [`Path`](@ref).
	"""

	width = adjuster.shift_width
	shift = rand(path.n_dimensions) * width .* rand([-1,1],path.n_dimensions)

	# Attempt one Single move
	adjuster.attempt_counter += 1

	# Save configuration of paths. Return to this configuration if Metropolis rejects the new configuration.
	old_beads = copy(path.beads[:, particle, :])

    # Evaluate old action. 
	# The kinetic action is unchanged because the relative bead positions are preserved, so we only need to evaluate the potential action at each bead.
    old_action = sum(potential_action(path, bead, particle, potential, regime) for bead in 1:path.n_beads)
	
	# Displace every bead on each specified particle by the same random vector.
	for bead in 1:path.n_beads
		path.beads[bead, particle, :] += shift
	end

    new_action = sum(potential_action(path, bead, particle, potential, regime) for bead in 1:path.n_beads)

	if new_action - old_action <= 0.0
		adjuster.success_counter += 1
		return true
	elseif rand() <= exp(-(new_action - old_action))
		adjuster.success_counter += 1
		return true
	else
		path.beads[:, particle, :] = old_beads
		return false
	end
end


function Bisect!(path::Path, particle::Int, potential::Potential, regime::Regime, adjuster::Adjuster)

	max_level = 3
	segment_length = Int((2^max_level) + 1)

	start_bead = rand(1:path.n_beads)
	old_beads = deepcopy(path.beads[:,particle, :])

	# Attempt one Single move
	adjuster.attempt_counter += 1

	total_old_action = 0.0
	for bead in start_bead:start_bead+segment_length
		total_old_action += potential_action(path, bead, particle, potential, regime)
	end

	for level in max_level:-1:1

		# Number of Beads in Level
		beads_level = 2^(max_level - level)

		for k in 1:beads_level

			# Find bead of interest
			if k == 1
				bead = Int(start_bead + (2^(level-1) * k))
			else
				bead = Int(start_bead + (2^(level-1) + (2^(level) * (k-1))))
			end

			# Find beads of which to find midpoint
			r0, r1 = bead - 2^(level-1), bead + 2^(level-1)

			#println("level: ", level, " bead: ", bead, " midpoints: ", [r0, r1], "\n")

			# Move by normally distributed shift
			#width = sqrt( 2^(level-1) * path.τ * path.λ) * adjuster.shift_width
			width = sqrt( 2^(level-1) * path.τ * path.λ) * adjuster.shift_width
			#shift = rand([-1,1],path.n_dimensions) .* rand(path.n_dimensions) * width
			shift = width .* rand(Distributions.Normal(0, 1), path.n_dimensions) # We want normal distribution
			
			# Perform move and calculate change to action
			midpoint = 0.5 * (path.beads[r0, particle, :] + path.beads[r1, particle, :])
			path.beads[bead, particle, :] = midpoint + shift
		end
	end

	print("\n")

	total_new_action = 0.0
	for bead in start_bead:start_bead+segment_length
		total_new_action += potential_action(path, bead, particle, potential, regime)
	end

	

	if total_new_action - total_old_action < 0.0
		adjuster.success_counter += 1
		#println("Success1")
		#println(path.beads[:, particle, :])
		return true

	elseif rand() < exp(-(total_new_action - total_old_action))
		adjuster.success_counter += 1
		#println("Success2")
		#println(path.beads[:, particle, :])
		return true
		
	else
		path.beads[:, particle, :] = old_beads
		#println(path.beads[:, particle, :])
		return false
	end

end


#=
function Bisect!(path::Path, particle::Int, potential::Potential, regime::Regime, adjuster::Adjuster)

	#max_level = Int(floor(log(rand(1:path.n_beads)) / log(2)))
	segment_length = 8 + 1 #temporary arbitary choice
	max_level = Int(floor(log(segment_length)/log(2))) # = 4 in arbitary setting

	start_bead = rand(1:path.n_beads)
	old_beads = deepcopy(path.beads[:,particle, :])

	# Attempt one Single move
	adjuster.attempt_counter += 1

	total_old_action = 0.0
	for bead in start_bead:start_bead+segment_length
		total_old_action += potential_action(path, bead, particle, potential, regime)
	end

	for level in max_level:-1:1
		segment_old_action = 0.0 # old action of the cut out segment
		segment_new_action = 0.0 # new action of the cut out segment

		ratio = 2^(max_level - level) #how many divisions of level makes up full segment
		
		for k in 1:ratio
			bead = Int(start_bead + (2^(level-1) * k))
			print("level: ", level, " bead: ", bead, "\n")
			#println("bead = ", bead)
			segment_old_action += potential_action(path, bead, particle, potential, regime)
			width = adjuster.shift_width * (2^0.5)^(level-1)
			width = sqrt( 2^(level-1) * path.τ * path.λ) # adjuster.shift_width[level]
			shift = rand([-1,1],path.n_dimensions) .* rand(path.n_dimensions) * width
			path.beads[bead, particle, :] = 0.5 * (path.beads[bead - 2^(level-1), particle, :] + path.beads[bead + 2^(level-1), particle, :]) + shift
			segment_new_action += potential_action(path, bead, particle, potential, regime)
		end
		segment_action_diff = 2^(level-1)* path.τ * (segment_new_action - segment_old_action)
		#=
		if rand() > exp(-segment_action_diff)
			print("Move Rejected")
			return false
		end
		=#
	end

	total_new_action = 0.0
	for bead in start_bead:start_bead+segment_length
		total_new_action += potential_action(path, bead, particle, potential, regime)
	end

	if total_new_action - total_old_action < 0.0
		adjuster.success_counter += 1
		return true

	elseif rand() < exp(-(total_new_action - total_old_action))
		adjuster.success_counter += 1
		return true
		
	else
		path.beads[:, particle, :] = old_beads
		return false
	end
end
=#