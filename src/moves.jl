# moves.jl
using Distributions
using TimerOutputs

#const tmr = TimerOutput();

function moveBead!(
    mover::SingleMover,
    path::Path,
    particle::Int,
    potential::Potential,
    regime::Regime,
    action_arr::Vector{Float64},
    shift::Vector{Float64},
    store_diff::Vector{Float64},
    prop_Matrix::Array{Float64},
)

    """
    Move a single imaginary-time timeslice (bead) on a single particle, or subset of particles, per Monte-Carlo iteration.

    Arguments
    - `mover::Mover: Single Mover type for this method`
    - `path::Path`: collection of all particle imaginary-time paths.
    - `particle::Int`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
    - `potential::Union{Potential}`: list of potentials active in the system. Can just be a single potential.
    See also [`Path`](@ref).
    """

    # Randomly choose a bead from the particle
    bead = rand(1:path.n_beads)
    width = mover.adjusters[particle].value

    #@timeit tmr "shift" 
    for i = 1:path.n_dimensions
        shift[i] = rand() * width * rand([-1, 1])# Linear random displacement of bead.
    end

    # Attempt one Single move
    mover.adjusters[particle].attempt_counter += 1

    # We just need to look at the kinetic contribution from beads +- 1 unit from the selected bead
    # Just calculate the potential action change at the selected bead
    #@timeit tmr "old_action" 
    action_arr[1] =
        kineticAction(path, bead - 1, bead, particle, regime, potential, store_diff) +# Link bead-1 to bead
        kineticAction(path, bead, bead + 1, particle, regime, potential, store_diff) +# Link bead to bead+1
        potentialAction(path, bead, particle, potential, regime, store_diff, prop_Matrix)# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.

    path.beads[bead, particle, :] += shift

    #@timeit tmr "new_action" 
    action_arr[2] =
        kineticAction(path, bead - 1, bead, particle, regime, potential, store_diff) +# Link bead-1 to bead
        kineticAction(path, bead, bead + 1, particle, regime, potential, store_diff) +# Link bead to bead+1
        potentialAction(path, bead, particle, potential, regime, store_diff, prop_Matrix)# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.


    # Metropolis algorithm. 
    # Accept if bead displacement decreases the action, otherwise accept with probability exp(-ΔS).
    # new - old action

    #@timeit tmr "decision" 
    if action_arr[2] - action_arr[1] <= 0.0 ||
       rand() <= exp(-(action_arr[2] - action_arr[1]))
        #if rand() <= exp(-(action_arr[2] - action_arr[1]))

        # Updating counter for adjustment of shift width
        mover.adjusters[particle].success_counter += 1
        return true

    else

        # Since rejected so we revert the shift
        path.beads[mod1(bead, path.n_beads), particle, :] -= shift
        return false

    end

end

function moveBead!(
    mover::DisplaceMover,
    path::Path,
    particle::Int,
    potential::Potential,
    regime::Regime,
)

    """
    Move the entire imaginary-time timeslice (all the beads) on a single particle, or subset of particles, by the same vector per Monte-Carlo iteration.

    Arguments
    Move a single imaginary-time timeslice (bead) on a single particle, or subset of particles, per Monte-Carlo iteration.

    Arguments
    - `mover::Mover: Displace Mover type for this method`
    - `path::Path`: collection of all particle imaginary-time paths.
    - `particle::Int`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
    - `potential::Union{Potential}`: list of potentials active in the system. Can just be a single potential.
    See also [`Path`](@ref).
    """

    width = mover.adjusters[particle].value
    shift = rand(path.n_dimensions) * width .* rand([-1, 1], path.n_dimensions)

    # Attempt one Single move
    mover.adjusters[particle].attempt_counter += 1

    # Save configuration of paths. Return to this configuration if Metropolis rejects the new configuration.
    old_beads = copy(path.beads[:, particle, :])

    # Evaluate old action. 
    # The kinetic action is unchanged because the relative bead positions are preserved, so we only need to evaluate the potential action at each bead.
    old_action = sum(
        potentialAction(path, bead, particle, potential, regime) for bead = 1:path.n_beads
    )

    # Displace every bead on each specified particle by the same random vector.
    for bead = 1:path.n_beads
        path.beads[bead, particle, :] += shift
    end

    new_action = sum(
        potentialAction(path, bead, particle, potential, regime) for bead = 1:path.n_beads
    )

    # Metropolis algorithm. 
    # Accept if bead displacement decreases the action, otherwise accept with probability exp(-ΔS).
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

function moveBead!(
    mover::BisectMover,
    path::Path,
    particle::Int,
    potential::Potential,
    regime::Regime,
)

    """
    Move a segment of imaginary-time timeslice (bead) on a single particle, or subset of particles, per Monte-Carlo iteration.
    Only works for simple potential (NOT for Frohlich model)

    ------Algorithm------
    1. Select a bead-segment length of L=2n+1 from a single particle
    2. Keep the start and end bead fixed
    3. "Remove" all the beads between the mid-points and construct a new path from mid-point
    4a. If such move accepted -> Move the mid-point beads in the subsequent level
    4b. If such move rejected -> Return the original position and deem the bisection move failed
    5. Repeat the mid-points displacement procedures until all levels of beads in the segment is accepted

    Arguments
    - `mover::Mover: Single Mover type for this method`
    - `path::Path`: collection of all particle imaginary-time paths.
    - `particle::Int`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
    - `potential::Potential`: list of potentials active in the system. Can just be a single potential.
    See also [`Path`](@ref).
    """

    # Segment length to perform Bisection
    max_level = mover.adjusters[particle].value # Level of adjustor
    segment_length = Int((2^max_level) + 1)

    # Make a copy of the old configurations in case of a rejection
    start_bead = rand(1:path.n_beads)
    old_beads = deepcopy(path.beads[:, particle, :])

    # Call parameters once
    τλ = path.τ * path.λ
    n_beads = path.n_beads
    n_dim = path.n_dimensions

    # Attempt one Single move
    mover.adjusters[particle].attempt_counter += 1

    # Calculates the old action
    total_old_action = 0.0
    for bead = start_bead:start_bead+segment_length
        total_old_action += potentialAction(path, bead, particle, potential, regime)
    end

    # Displacement of mid-points in each level
    for level = max_level:-1:1
        level_fac = 2^(level - 1)

        # Number of Beads in Level
        beads_level = 2^(max_level - level)

        for k = 1:beads_level

            # Find bead of interest
            if k == 1
                bead = Int(start_bead + (level_fac * k))
            else
                bead = Int(start_bead + (level_fac + (2^(level) * (k - 1))))
            end

            # Find beads of which to calculates the midpoint coordinates
            r0, r1 = bead - level_fac, bead + level_fac

            # Find normally distributed shift
            shift = sqrt(level_fac * τλ) .* rand(Distributions.Normal(0, 1), n_dim)

            # Calculate the midpoint & perform move
            path.beads[mod1(bead, n_beads), particle, :] =
                0.5 * (
                    path.beads[mod1(r0, n_beads), particle, :] +
                    path.beads[mod1(r1, n_beads), particle, :]
                ) + shift
        end
    end

    # Calculates the new action (and ΔS)
    total_new_action = 0.0
    for bead = start_bead:start_bead+segment_length
        total_new_action += potentialAction(path, bead, particle, potential, regime)
    end

    # Determine whether to accept the full move
    if total_new_action - total_old_action <= 0.0 ||
       rand() <= exp(-(total_new_action - total_old_action))
        mover.adjusters[particle].success_counter += 1

        return true

    else
        path.beads[:, particle, :] = old_beads

        return false
    end
end

function moveBead!(
    mover::BisectMover,
    path::Path,
    particle::Int,
    potential::FrohlichPotential,
    regime::Regime,
)
    """
    Special Version of Bisection algorithm for Frohlich Potential
    as Frohlich potential contains problem with double counting contributions when calculating potential changes

    """
    # Segment length to perform Bisect
    max_level = mover.adjusters[particle].value
    segment_length = Int((2^max_level))

    # Randomly choose a starting bead
    start_bead = rand(1:path.n_beads)
    old_beads = deepcopy(path.beads[:, particle, :])

    # Call parameters once to save time when calling struct attributes
    τλ = path.τ * path.λ
    n_beads = path.n_beads
    n_dim = path.n_dimensions

    # Define the beadrange to avoid double counting beads between different layer
    beadrange = collect(start_bead+1:start_bead+segment_length-1)
    for i = 1:segment_length-1
        beadrange[i] = mod1(beadrange[i], n_beads)
    end

    # Attempt one Bisect move. Only accept when the whole segment get accepted, otherwise count as rejection
    mover.adjusters[particle].attempt_counter += 1

    # Calculates the total old_action
    total_old_action = 0.0
    for bead = start_bead+1:start_bead+segment_length-1
        total_old_action +=
            potentialAction(path, bead, beadrange, particle, potential, regime)
    end

    # Displace the beads level by level respectively. Only proceed to the next layer if the previous layer is accepted
    for level = max_level:-1:1
        level_fac = 2^(level - 1)

        # Number of Beads in Level
        beads_level = 2^(max_level - level)

        for k = 1:beads_level

            # identify bead of interest
            if k == 1
                bead = Int(start_bead + (level_fac * k))
            else
                bead = Int(start_bead + (level_fac + (2^(level) * (k - 1))))
            end

            # Find beads of which to find midpoint
            r0, r1 = bead - level_fac, bead + level_fac

            # Shifting the beads based on 1/2*[r0+r1] + sqrt(level * τλ) * normal_distribution(0, 1) [Ref TESI paper]
            path.beads[mod1(bead, n_beads), particle, :] =
                0.5 * (
                    path.beads[mod1(r0, n_beads), particle, :] +
                    path.beads[mod1(r1, n_beads), particle, :]
                ) + sqrt(level_fac * τλ) .* rand(Distributions.Normal(0, 1), n_dim)

            # Infinite potential well (Can put a artifically large well_size if not required)
            if any(abs.(path.beads[mod1(bead, n_beads), particle, :]) .> well_size)
                path.beads[:, particle, :] = old_beads
                return false
            end
        end
    end

    # Calculates new action
    total_new_action = 0.0
    for bead = start_bead+1:start_bead+segment_length-1
        total_new_action +=
            potentialAction(path, bead, beadrange, particle, potential, regime)
    end


    # Metropolis algorithm, decide accept or reject
    if total_new_action - total_old_action <= 0.0 ||
       rand() <= exp(-(total_new_action - total_old_action))
        mover.adjusters[particle].success_counter += 1

        return true

    else
        path.beads[:, particle, :] = old_beads

        return false
    end
end

function moveBead!(
    path::DiscretePath,
    particle::Int,
    potential::Potential,
    F_l::Array{Float64},
)

    """
    Discretised single bead movement
    Move a single imaginary-time timeslice (bead) on a single particle, or subset of particles, per Monte-Carlo iteration.

    Arguments
    - `mover::Mover: Single Mover type for this method`
    - `path::Path`: collection of all particle imaginary-time paths.
    - `particle::Int`: select a specific particle indexed by this integer, or a subset of particles indexed by integers in this vector.
    - `potential::Union{Potential}`: list of potentials active in the system. Can just be a single potential.
    See also [`Path`](@ref).
    """

    # Randomly choose a bead from the particle
    bead = rand(1:path.n_beads)
    dimension = rand(1:path.n_dimensions)
    shift = rand([-1, 1]) # Random shifting

    # We just need to look at the kinetic contribution from beads +- 1 unit from the selected bead
    # Just calculate the potential action change at the selected bead
    #@timeit tmr "old_action" 
    old_action =
        kineticAction(path, bead, mod1(bead - 1, path.n_beads), particle, potential) *
        kineticAction(path, bead, mod1(bead + 1, path.n_beads), particle, potential) *
        exp(potentialAction(path, bead, particle, potential, F_l))# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.

    path.beads[bead, particle, dimension] += shift
    #path.beads[bead] += shift

    #@timeit tmr "new_action" 
    new_action =
    #prod([kineticAction(path, bead, other_bead, particle, potential) for other_bead in 1:path.n_beads if other_bead != bead]) * 
        kineticAction(path, bead, mod1(bead - 1, path.n_beads), particle, potential) *
        kineticAction(path, bead, mod1(bead + 1, path.n_beads), particle, potential) *
        exp(potentialAction(path, bead, particle, potential, F_l))# Potential at bead for all particles incl. any const., 1-body or 2-body interactions.


    # Metropolis algorithm. 
    # Accept if bead displacement decreases the action, otherwise accept with probability exp(-ΔS).
    # new - old action
    #println("NEW:", trunc(new_action, digits=2), " ", "OLD:", trunc(old_action, digits=2))
    #println("ratio:", trunc(new_action/old_action, digits=2))

    if rand() <= new_action / old_action
        # Updating counter for adjustment of shift width
        #mover.adjusters[particle].success_counter += 1
        return true

    else
        # Since rejected so we revert the shift
        path.beads[mod1(bead, path.n_beads), particle, dimension] -= shift
        #path.beads[mod1(bead, path.n_beads)] -= shift
        return false

    end

end
