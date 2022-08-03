# moves.jl




#Sampling method of moving single beads individually
function Single!(path::Path, particle::Int, potential::Potential, regime::Regime, adjuster::Adjuster)
    bead = rand(1:path.n_beads)
	width = adjuster.shift_width
	shift = width * rand([-1,1]) * rand()

    # We just need to look at the beads +- 1 unit from m
    # CHECK: Non local potential? Coulombic?
    old_action = 
		kinetic_action(path, bead-1, bead, particle, regime) +		# Link bead-1 to bead
		kinetic_action(path, bead, bead+1, particle, regime) +		# Link bead to bead+1
		potential_action(path, bead, particle, potential, regime)	# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.

	path.beads[bead, particle] += shift

    new_action =
		kinetic_action(path, bead-1, bead, particle, regime) +		# Link bead-1 to bead
		kinetic_action(path, bead, bead+1, particle, regime) +		# Link bead to bead+1
		potential_action(path, bead, particle, potential, regime)	# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.


	if new_action - old_action <= 0.0 || rand() <= exp(-(new_action - old_action))
		adjuster.adjust_counter += 1 #updating counter for adjustment of shift width
		return true
	else
		path.beads[bead, particle] -= shift
		adjuster.adjust_counter -= 1 #updating counter for adjustment of shift width
		return false
	end
end




function Displace!(path::Path, particle::Int, potential::Potential, regime::Regime, adjuster::Adjuster)
	width = adjuster.shift_width
	shift = width * rand([-1,1]) * rand()

	old_beads = copy(path.beads[:, particle])
    
    old_action = sum(potential_action(path, bead, particle, potential, regime) for bead in 1:path.n_beads)
	
	for bead in 1:path.n_beads
		path.beads[bead, particle] += shift
	end

    new_action = sum(potential_action(path, bead, particle, potential, regime) for bead in 1:path.n_beads)

	if new_action - old_action <= 0.0
		adjuster.adjust_counter += 1 #updating counter for adjustment of shift width
		return true
	elseif rand() <= exp(-(new_action - old_action))
		adjuster.adjust_counter += 1 #updating counter for adjustment of shift width
		return true
	else
		path.beads[:, particle] = old_beads
		adjuster.adjust_counter -= 1 #updating counter for adjustment of shift width
		return false
	end
end


function Bisect!(path::Path, particle::Int, potential::Potential, regime::Regime, adjuster::Adjuster)

	#max_level = Int(floor(log(rand(1:path.n_beads)) / log(2)))
	segment_length = 16 + 1 #temporary arbitary choice
	max_level = Int(floor(log(segment_length)/log(2))) # = 4 in arbitary setting

	start_bead = rand(1:path.n_beads)

	old_beads = copy(path.beads[:,particle])


	total_old_action = 0.0
	for bead in start_bead:start_bead+segment_length-1
		total_old_action += total_action(path, bead, bead+1, particle, potential, regime)
	end


	segment_old_action = 0.0 # old action of the cut out segment
	segment_new_action = 0.0 # new action of the cut out segment
	for level in max_level-1:-1:0

		
		ratio = (segment_length - 1) / 2^level #how many divisions of level makes up full segment

		for interval in 1:ratio-1
			bead = Int(start_bead + (2^level * interval))
			segment_old_action += total_action(path, bead, bead+1, particle, potential, regime)
			shift = adjuster.shift_width_array[string(level)]*2 #auto adjusting level specific width
			path.beads[bead, particle] = 0.5 * (path.beads[bead - 2^level, particle] + path.beads[bead + 2^level, particle]) + shift
			segment_new_action += total_action(path, bead, bead+1, particle, potential, regime)
		end
		segment_action_diff = segment_new_action - segment_old_action
		if rand() >= exp(-segment_action_diff)
			adjuster.adjust_counter_array[string(level)] -= 1
			return false
		else
			adjuster.adjust_counter_array[string(level)] += 1

		end

	end

	total_new_action = 0.0
	for bead in start_bead:start_bead+segment_length-1
		total_new_action += total_action(path, bead, bead+1, particle, potential, regime)
	end

	if total_new_action - total_old_action < 0.0

		return true
	elseif rand() < exp(-(total_new_action - total_old_action))

		return true
	else
		path.beads[:, particle] = old_beads
		return false
	end


end







