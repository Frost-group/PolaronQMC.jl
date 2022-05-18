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
    β = path.τ * path.n_beads
    inner_integral = 0.0
    for other_bead in 1:path.n_beads
        if other_bead != bead
            g_factor = potential.α/2 * sqrt(potential.ħ/(2*path.m*potential.ω)) * cosh(β * (abs(bead-other_bead)/path.n_beads - 0.5)) * csch(β/2)
            inner_integral += g_factor / norm(path.beads[bead, particle, :] - path.beads[other_bead, particle, :])
        end
    end
    return path.τ * inner_integral
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

abstract type TwoBodyPotential <: Potential end

struct ZeroPotential <: Potential end

struct HarmonicPotential <: OneBodyPotential
    ω :: Float64
    function HarmonicPotential(ω::Float64)
        new(ω)
    end
end

struct CoulombPotential <: TwoBodyPotential
    κ :: Float64
    function CoulombPotential(κ::Float64)
        new(κ)
    end
end

function one_body_potential(potential::ZeroPotential, path::Path, bead::Int, particle::Int)
    return 0.0
end

function one_body_potential(potential::HarmonicPotential, path::Path, bead::Int, particle::Int)
    return 0.5 * potential.ω^2 * norm(path.beads[bead, particle, :])^2
end

function two_body_potential(potential::ZeroPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
    return 0.0
end

function two_body_potential(potential::CoulombPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
    return -potential.κ / norm(path.beads[bead, particle_one, :] .- path.beads[bead, particle_two, :])
end
