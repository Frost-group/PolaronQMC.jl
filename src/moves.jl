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
function Single!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})
    bead = rand(1:path.n_beads)								# Pick a random bead.
	width = sqrt(2 * path.λ * path.τ)						# Displacement width. ~Order(thermal de Broglie wavelength). Adjust for ~50% acceptance rate.
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
	if new_action - old_action <= 0.0 || rand() <= exp(-(new_action - old_action))
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
function Displace!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})
	width = sqrt(150 * path.λ * path.τ)			# Displacement width. ~O(thermal de Broglie wavelength). Adjust for ~50% acceptance rate.
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
	if new_action - old_action <= 0.0 || rand() <= exp(-(new_action - old_action))
		return true
	else
		path.beads[:, particle, :] = old_beads	# New path configuration rejected, return to old configuration.
		return false
	end
end

function Stage!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})
	relabel_beads!(path)

	segment_length = 16

	old_action = 0.0
	for bead in 2:segment_length-1
		old_action += potential_action(path, bead, particle, potentials)
	end

	old_beads = copy(path.beads[:, particle, :])
	new_action = 0.0

	for bead in 1:segment_length-1
		staging_mass = (segment_length - bead + 1) / (segment_length - bead)
		staging_position = (path.beads[1 + segment_length, particle, :] + path.beads[bead, particle, :] * (segment_length - bead)) / (segment_length - bead + 1)
		path.beads[bead + 1, particle, :] += staging_position + randn(path.n_dimensions) * sqrt(path.τ / staging_mass)
		new_action += potential_action(path, bead + 1, particle, potentials)
	end

	if rand() <= exp(-(new_action - old_action))
		return true
	else
		path.beads[:, particle, :] = old_beads
		return false
	end
end

function Bisect!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})
	relabel_beads!(path)

	max_level = Int(floor(log(rand(1:path.n_beads)) / log(2)))
	clip_length = 2^max_level + 1

	old_action = 0.0
	for bead in 1:clip_length
		old_action += primitive_action(path, bead, particle, potentials)
	end

	for level in max_level:-1:1
		step = 2^(level - 1)
		shift = randn(path.n_dimensions) * sqrt(step * path.τ * path.λ)

		old_action = 
		
		for n in 1:2^(max_level - level)
			bead = 1 + n * step
			old_action = primitive_action(path, bead, particle, potentials)
			path.beads[bead, particle, :] = (path.beads[bead - step, particle, :] + path.beads[bead + step, particle, :]) / 2 + shift
		end
	end

	old_action = primitive_action(path, 1, particle, potentials)
end

