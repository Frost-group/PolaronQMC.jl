# potentials.jl

"""
Collection of potentials to use in simulations. 
"""

"""
Outer constructors for different potential types.
"""

# Just return value of potential for a constant potential independent of single particle.
function one_body_potential(potential::ConstantPotential, potentialcache::Cache, path::Path, bead::Int, particle::Int)
    return potential.V
end

# Return the harmonic potential for a single particle.
function one_body_potential(potential::HarmonicPotential, potentialcache::Cache, path::Path, bead::Int, particle::Int)
    return potentialcache.cached_value * norm(path.beads[bead, particle,:])^2
end

# Returns the Frohlich potential for a single particle
function one_body_potential(potential::FrohlichPotential, potentialcache::Cache, path::Path, bead::Int, particle::Int)
    β = path.n_beads * path.τ
    bead_distances = filter!(!iszero, potentialcache.distance_matrix[bead,:])
    
    inner_integral = sum(potentialcache.potential_prefactor[bead] ./ bead_distances)
    return path.τ * inner_integral
end


function one_body_potentialL(potential::FrohlichPotential, potentialcache::Cache, path::Path, bead::Int, particle::Int)
    β = path.n_beads * path.τ
    inner_integral = 0.0
    for other_bead in 1:path.n_beads
        if other_bead != bead
            #g_factor = -potential.α/2 * sqrt(potential.ħ/(2*path.m*potential.ω)) * cosh(β *(abs(bead-other_bead)/path.n_beads - 0.5)) * csch(β/2)
            g_factor = -0.5 * potential.α * (potential.ħ * potential.ω)^3/2 * sqrt(2*path.m) * csch(potential.ħ * potential.ω * β / 2) * cosh(potential.ω*β * (abs(bead-other_bead)/path.n_beads - 0.5 * potential.ħ)) 
            inner_integral += g_factor / norm(path.beads[bead, particle, :] - path.beads[other_bead, particle, :])
        end
    end
    return path.τ * inner_integral
end


# Mexican Hat -r^2+r^4 in N-dimensions
function one_body_potential(potential::MexicanHatPotential, potentialcache::Cache, path::Path, bead::Int, particle::Int)
    r=norm(path.beads[bead,particle,:])^2
    return 0.5 * potential.ω^2 * (-r^2+r^4)
end

# Just return value of potential for a constant potential independent of two particles.
function two_body_potential(potential::ConstantPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
    return potential.V
end

# Return the Coulomb potential between two particles.
function two_body_potential(potential::CoulombPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
    return -potential.κ / norm(path.beads[bead, particle_one, :] .- path.beads[bead, particle_two, :])
end

