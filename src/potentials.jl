# potentials.jl

"""
Collection of potentials to use in simulations. 
"""

"""
Outer constructors for different potential types.
"""

# Just return value of potential for a constant potential independent of single particle.
function one_body_potential(potential::ConstantPotential, path::Path, bead::Int, particle::Int)
    return potential.V
end

# Return the harmonic potential for a single particle.
function one_body_potential(potential::HarmonicPotential, path::Path, bead::Int, particle::Int)
    return 0.5 * path.m * potential.ω^2 * norm(path.beads[bead, particle,:])^2
end

# Returns the Frohlich potential for a single particle
function one_body_potential(potential::FrohlichPotential, path::Path, bead::Int, particle::Int)
    ω = potential.ω
    phonon_response(other_bead) = cosh(ω * path.τ * (abs(bead - other_bead) - path.n_beads / 2.0)) / sinh(ω * path.τ * path.n_beads / 2.0)
    double_integral = 0.0
    for other_bead in 1:path.n_beads
        if other_bead != bead
            double_integral += phonon_response(other_bead) / norm(path.beads[bead, particle, :] .- path.beads[other_bead, particle, :])
        end
    end
    return potential.α * ω^(-1/2) / sqrt(8) * double_integral / path.n_beads
end



# Mexican Hat -r^2+r^4 in N-dimensions
function one_body_potential(potential::MexicanHatPotential, path::Path, bead::Int, particle::Int)
    r=norm(path.beads[bead,particle])^2
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

