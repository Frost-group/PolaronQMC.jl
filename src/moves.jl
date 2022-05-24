# moves.jl

function Single!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})
    bead = rand(1:path.n_beads)
	width = sqrt(path.λ * path.τ)
	shift = width .* (2 .* rand(path.n_dimensions) .- 1)

    # We just need to look at the beads +- 1 unit from m
    # CHECK: Non local potentials? Coulombic?
    old_action = 
		kinetic_action(path, bead-1, bead, particle) +		# Link bead-1 to bead
		kinetic_action(path, bead, bead+1, particle) +		# Link bead to bead+1
		potential_action(path, bead, particle, potentials)	# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.

	path.beads[bead, particle, :] += shift

    new_action =
		kinetic_action(path, bead-1, bead, particle) +		# Link bead-1 to bead
		kinetic_action(path, bead, bead+1, particle) +		# Link bead to bead+1
		potential_action(path, bead, particle, potentials)	# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.


	if new_action - old_action <= 0.0 || rand() <= exp(-(new_action - old_action))
		return true
	else
		path.beads[bead, particle, :] -= shift
		return false
	end
end

function Displace!(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})
	width = sqrt(path.λ * path.τ)
	shift = width .* randn(path.n_dimensions)

	old_beads = copy(path.beads[:, particle, :])
    
    old_action = sum(potential_action(path, bead, particle, potentials) for bead in 1:path.n_beads)
	
	for bead in 1:path.n_beads
		path.beads[bead, particle, :] += shift
	end

    new_action = sum(potential_action(path, bead, particle, potentials) for bead in 1:path.n_beads)

	if new_action - old_action <= 0.0
		return true
	elseif rand() <= exp(-(new_action - old_action))
		return true
	else
		path.beads[:, particle, :] = old_beads
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

