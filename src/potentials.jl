# potentials.jl

abstract type Potential end

abstract type OneBodyPotential <: Potential end

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
