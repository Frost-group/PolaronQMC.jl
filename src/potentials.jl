# potentials.jl

"""
Collection of potentials to use in simulations. 
"""

"""
Outer constructors for different potential types.
"""

# Just return value of potential for a constant potential independent of single particle.
function oneBodyPotential(potential::ConstantPotential, path::Path, bead::Int, particle::Int)
    return potential.V
end

# Return the harmonic potential for a single particle.
function oneBodyPotential(potential::HarmonicPotential, path::Path, bead::Int, particle::Int)
    return 0.5 * path.m * potential.ω^2 * norm(path.beads[mod1(bead, path.n_beads), particle,:])^2
end


function oneBodyPotential(potential::HarmonicInteractionPotential, path::Path, bead::Int, particle::Int)
    harmonic_pot = 0.5 * path.m * potential.ω^2 * norm(path.beads[mod1(bead, path.n_beads), particle,:])^2
    coloumb_int = 0
    for other_particle in 1:path.n_particles
        if other_particle != particle
            coloumb_int += -potential.κ / norm(path.beads[mod1(bead, path.n_beads), particle, :] .- path.beads[mod1(bead, path.n_beads), other_particle, :])
        end
    end
    return harmonic_pot + coloumb_int
end

# Returns the Frohlich potential for a single particle

function oneBodyPotential(potential::FrohlichPotential, path::Path, bead::Int, particle::Int)
    β = path.τ * path.n_beads
    m = path.m
    ħω = potential.ħ * potential.ω
    α = potential.α
    term_factor = -0.5 * α * (ħω)^(3/2) * sqrt(1/2/m) * csch(ħω * β / 2)
    
    inner_integral = 0.0
    for other_bead in 1:path.n_beads
        if other_bead != bead
            inner_integral += cosh(ħω * β * (abs(bead-other_bead)/path.n_beads - 0.5)) / norm(path.beads[mod1(bead, path.n_beads), particle, :] - path.beads[mod1(other_bead, path.n_beads), particle, :])
        end
    end
    return path.τ * inner_integral * term_factor # Note that this path.τ multiplication refer to dτ'
end


function oneBodyPotential(potential::FrohlichInteractionPotential, path::Path, bead::Int, particle::Int)
    β = path.τ * path.n_beads
    m = path.m
    ħω = potential.ħ * potential.ω
    α = potential.α
    term_factor = -0.5 * α * (ħω)^(3/2) * sqrt(1/2/m) * csch(ħω * β / 2)
    
    inner_integral = 0.0
    for other_bead in 1:path.n_beads
        if other_bead != bead
            inner_integral += cosh(ħω * β * (abs(bead-other_bead)/path.n_beads - 0.5)) / norm(path.beads[mod1(bead, path.n_beads), particle, :] - path.beads[mod1(other_bead, path.n_beads), particle, :])
        end
    end

    coloumb_int = 0
    for other_particle in 1:path.n_particles
        if other_particle != particle
            coloumb_int += -potential.κ / norm(path.beads[mod1(bead, path.n_beads), particle, :] .- path.beads[mod1(bead, path.n_beads), other_particle, :])
        end
    end
    return path.τ * inner_integral * term_factor + coloumb_int # Note that this path.τ multiplication refer to dτ'
end


# Mexican Hat -r^2+r^4 in N-dimensions
function oneBodyPotential(potential::MexicanHatPotential, path::Path, bead::Int, particle::Int)
    r = norm(path.beads[mod1(bead, path.n_beads), particle,:])^2
    return 0.5 * potential.ω^2 * (-r^2+r^4)
end


# Just return value of potential for a constant potential independent of two particles.
function twoBodyPotential(potential::ConstantPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
    return potential.V
end


# Return the Coulomb potential between two particles.
function twoBodyPotential(potential::CoulombPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
    return -potential.κ / norm(path.beads[mod1(bead, path.n_beads), particle_one, :] .- path.beads[mod1(bead, path.n_beads), particle_two, :])
end

