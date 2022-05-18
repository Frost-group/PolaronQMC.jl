# potentials.jl

"""
Collection of potentials to use in simulations. 
"""

"""
Defined Types for different potentials.
"""

# Most abstract Potential type.
abstract type Potential end

# A potential that is independent on bodies.
abstract type NoBodyPotential <: Potential end

# For potentials that depend on one body.
abstract type OneBodyPotential <: Potential end

# For potentials that depend on two bodies.
abstract type TwoBodyPotential <: Potential end

# A constant potential.
struct ConstantPotential <: NoBodyPotential 
    V :: Float64
    function ConstantPotential(V::Float64)
        new(V)
    end
end

# A Harmonic potential for a single body.
struct HarmonicPotential <: OneBodyPotential
    ω :: Float64
    function HarmonicPotential(ω::Float64)
        new(ω)
    end
end

# Coulomb interaction between two bodies.
struct CoulombPotential <: TwoBodyPotential
    κ :: Float64
    function CoulombPotential(κ::Float64)
        new(κ)
    end
end

"""
Outer constructors for different potential types.
"""

# Just return value of potential for a constant potential independent of single particle.
function one_body_potential(potential::ConstantPotential, path::Path, bead::Int, particle::Int)
    return potential.V
end

# Return the harmonic potential for a single particle.
function one_body_potential(potential::HarmonicPotential, path::Path, bead::Int, particle::Int)
    return 0.5 * potential.ω^2 * norm(path.beads[bead, particle, :])^2
end

# Just return value of potential for a constant potential independent of two particles.
function two_body_potential(potential::ConstantPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
    return potential.V
end

# Return the Coulomb potential between two particles.
function two_body_potential(potential::CoulombPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
    return -potential.κ / norm(path.beads[bead, particle_one, :] .- path.beads[bead, particle_two, :])
end
