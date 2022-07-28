# moves.jl

"""
	Single!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})

Move a single imaginary-time timeslice (bead) on a single particle, or subset of particles, per Monte-Carlo iteration.

# Arguments
- `path::Path`: collection of all particle imaginary-time paths.
- `particle::Union{Int, Vector{Int}}`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
- `potentials::Union{Potential, Array{Potential}}`: list of potentials active in the system. Can just be a single potential.

See also [`Path`](@ref).
"""
function Single!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}}; scale = 4.0)
    bead = rand(1:path.n_beads)								# Pick a random bead.
	width = sqrt(scale * path.λ * path.τ)						# Displacement width. ~Order(thermal de Broglie wavelength). Adjust for ~50% acceptance rate.
	shift = width .* (2 .* rand(path.n_dimensions) .- 1)	# Linear random displacement of bead.	

	# Save configuration of paths. Return to this configuration if Metropolis rejects the new configuration.
	old_beads = copy(path.beads[bead, particle, :])

    # Evaluate action for links connected to selected bead for specified particle. I.e bead-1 to bead and bead to bead+1.
    old_action = 														
		primitive_action(path, bead-1, bead, particle, potentials) +	# Link bead-1 to bead.
		primitive_action(path, bead, bead+1, particle, potentials)		# Link bead to bead+1.

	# Displace the bead for specified particle.
	path.beads[bead, particle, :] += shift

	# Evaluate new action.
    new_action =														
		primitive_action(path, bead-1, bead, particle, potentials) +	# Link bead-1 to bead.
		primitive_action(path, bead, bead+1, particle, potentials)		# Link bead to bead+1.

	# Metropolis algorithm. 
	# Accept if bead displacement decreases the action, otherwise accept with probability exp(-ΔAction).
	if rand() < minimum([1,exp(-(new_action - old_action))])
		return true
	else
		path.beads[bead, particle, :] = old_beads	# Displacement rejected so return bead to prior position.
		return false
	end
end

"""
	Displace!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})

Move the entire imaginary-time timeslice (all the beads) on a single particle, or subset of particles, by the same vector per Monte-Carlo iteration.

# Arguments
- `path::Path`: collection of all particle imaginary-time paths.
- `particle::Union{Int, Vector{Int}}`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
- `potentials::Union{Potential, Array{Potential}}`: list of potentials active in the system. Can just be a single potential.

See also [`Path`](@ref).
"""
function Displace!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}}; scale = 4.0)
	width = sqrt(scale * path.λ * path.τ)			# Displacement width. ~O(thermal de Broglie wavelength). Adjust for ~50% acceptance rate.
	shift = width .* randn(path.n_dimensions)	# Normally distributed random displacement of particle.

	# Save configuration of paths. Return to this configuration if Metropolis rejects the new configuration.
	old_beads = copy(path.beads[:, particle, :])
    
	# Evaluate old action. 
	# The kinetic action is unchanged because the relative bead positions are preserved, so we only need to evaluate the potential action at each bead.
    old_action = sum(potential_action(path, bead, particle, potentials) for bead in 1:path.n_beads)
	
	# Displace every bead on each specified particle by the same random vector.
	for bead in 1:path.n_beads
		path.beads[bead, particle, :] += shift
	end

	# Evaluate new action.
    new_action = sum(potential_action(path, bead, particle, potentials) for bead in 1:path.n_beads)

	# Metropolis algorithm.
	# Accept if bead displacements decreases the potential action, otherwise accept with probability exp(-ΔPotentialAction).
	if rand() < minimum([1,exp(-(new_action - old_action))])
		return true
	else
		path.beads[:, particle, :] = old_beads	# New path configuration rejected, return to old configuration.
		return false
	end
end


function Bisect!(path::Path, particle::Int, potential::Potential, regime::Regime, adjuster::Adjuster)

	#max_level = Int(floor(log(rand(1:path.n_beads)) / log(2)))
	segment_length = 16 + 1 #temporary arbitary choice
	max_level = Int(floor(log(segment_length)/log(2))) # = 4 in arbitary setting

	start_bead = rand(1:path.n_beads)
	#println("sb = ", start_bead)
	#println("ed = ",start_bead+segment_length)

	old_beads = copy(path.beads[:,particle, :])


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
			#println("bead = ", bead)
			segment_old_action += potential_action(path, bead, particle, potential, regime)
			shift = rand([-1,1]) * rand(path.n_dimensions) * sqrt( 2^(level-1) * path.τ * path.λ) 
			path.beads[bead, particle, :] = 0.5 * (path.beads[bead - 2^(level-1), particle, :] + path.beads[bead + 2^(level-1), particle, :]) + shift
			segment_new_action += potential_action(path, bead, particle, potential, regime)
		end
		segment_action_diff = 2^(level-1)* path.τ * (segment_new_action - segment_old_action)
		if rand() > exp(-segment_action_diff)
			#adjuster.adjust_counter_array[string(level)] -= 1
			return false
		end

	end

	total_new_action = 0.0
	for bead in start_bead:start_bead+segment_length
		total_new_action += potential_action(path, bead, particle, potential, regime)
	end

	if total_new_action - total_old_action < 0.0

		return true
	elseif rand() < exp(-(total_new_action - total_old_action))

		return true
	else
		path.beads[:, particle, :] = old_beads
		return false
	end

end






