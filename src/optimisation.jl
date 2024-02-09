# optimisation.jl
"""
Algorithm for speeding up the simulation
1. Adjusting the shift width of beads by calculating the acceptance ratio
2. Adjusting the maximum number of beads shifted in bisect mover algorithm

"""


function updateAdjuster(adjuster::Union{SingleAdjuster,DisplaceAdjuster}, path::Path)
    """
    Update shift width by keeping track the acceptance ratio = successful shift/total shift
        (total shift = number of moves in one sweep) -> reset the ratio when a new sweep initiated
    """
    if (adjuster.attempt_counter != 0)
        adjuster.acceptance_rate = adjuster.success_counter / adjuster.attempt_counter

        if adjuster.acceptance_rate < 0.01 # To quickly lower the shift width if the acceptance rate is too low
            adjuster.value *= 0.99

        elseif adjuster.acceptance_rate > 0.99 # To quickly increase the shift width if the acceptance rate is too low
            adjuster.value *= 1.01

        else
            adjuster.value *= (adjuster.acceptance_rate / 0.3) # originally it's 0.5
        end

        # Resetting the count for next sweep
        adjuster.attempt_counter = 0
        adjuster.success_counter = 0
    end
end


function updateAdjuster(adjuster::Union{BisectAdjuster}, path::Path)
    """
    Update the bisect segment length based on acceptance rate
        if acceptance rate too low -> Decrease the number of bisecting level
        if acceptance rate too high -> Increase the number of bisecting level so can displace more beads at once
    """
    adjuster.acceptance_rate = adjuster.success_counter / adjuster.attempt_counter

    if adjuster.acceptance_rate < 0.5
        if adjuster.value > 0
            adjuster.value -= 1
        end

    else
        if 2^(adjuster.value + 1) + 1 < path.n_beads
            adjuster.value += 1
        end

    end

    # Resetting the count for next sweep
    adjuster.attempt_counter = 0
    adjuster.success_counter = 0
end

function copyLastPath!(
    path::Path,
    potential::Potential,
    A::Union{SizedArray,Array};
    verbose::Bool = true,
)
    path.beads[:, :, :] = A[:, :, :]
    if verbose
        println("Copying last path complete")
    end
end

function recentralise(path::Path; verbose = true)
    """
    Shifting the bead of chain for each particle such that their average position is 0
    To avoid polaron moving to infinities
    """
    for particle = 1:path.n_particles
        centroid_pos =
            [
                sum(path.beads[bead, particle, dimension] for bead = 1:path.n_beads) for
                dimension = 1:path.n_dimensions
            ] / path.n_beads
        for bead = 1:path.n_beads
            path.beads[bead, particle, :] -= centroid_pos
        end
    end
    if verbose
        println("Readjusted centre")
    end
end
