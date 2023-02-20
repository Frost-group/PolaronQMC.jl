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
    return 0.5 * path.m * potential.ω^2 * norm(path.beads[mod1(bead, path.n_beads), particle,:])^2
end

# Returns the Frohlich potential for a single particle
function one_body_potential(potential::FrohlichPotential, path::Path, bead::Int, particle::Int)
    β = path.τ * path.n_beads
    m = path.m
    ħω = potential.ħ * potential.ω
    α = potential.α
    term_factor = -0.5 * α * (ħω)^(3/2) * sqrt(1/2/m) * csch(ħω * β / 2)
    
    inner_integral = 0.0
    for other_bead in 1:path.n_beads
        if other_bead != bead
            #g_factor = -potential.α/2 * sqrt(potential.ħ/(2*path.m*potential.ω)) * cosh(β *(abs(bead-other_bead)/path.n_beads - 0.5)) * csch(β/2)
            #g_factor = -0.5 * potential.α * (potential.ħ * potential.ω)^3/2 * sqrt(2*path.m) * cosh(potential.ω*β * (abs(bead-other_bead)/path.n_beads - 0.5 * potential.ħ)) * csch(potential.ħ * potential.ω * β / 2)

            #g_factor = -0.5 * potential.α * (potential.ħ * potential.ω)^3/2 / pi * sqrt(2/path.m) * cosh(potential.ω * β * potential.ħ * (abs(bead-other_bead)/path.n_beads - 0.5))* csch(potential.ħ * potential.ω * β / 2)
            #g_factor = -0.5 * potential.α * (potential.ħ * potential.ω)^3/2 * sqrt(1/2/path.m) * cosh(potential.ω * β * potential.ħ * (abs(bead-other_bead)/path.n_beads - 0.5))* csch(potential.ħ * potential.ω * β / 2)
            #inner_integral += g_factor / norm(path.beads[bead, particle, :] - path.beads[other_bead, particle, :])
            #g_factor = -0.5 * α * (ħω)^(3/2) * sqrt(1/2/m) * cosh(ħω * β * (abs(bead-other_bead)/path.n_beads - 0.5))* csch(ħω * β / 2)
            #g_factor = cosh(ħω * β * (abs(bead-other_bead)/path.n_beads - 0.5))
            
            nner_integral += cosh(ħω * β * (abs(bead-other_bead)/path.n_beads - 0.5)) / norm(path.beads[mod1(bead, path.n_beads), particle, :] - path.beads[mod1(other_bead, path.n_beads), particle, :])
            #inner_integral += g_factor / norm(path.beads[mod1(bead, path.n_beads), particle, :] - path.beads[mod1(other_bead, path.n_beads), particle, :])
        end
    end
    return path.τ * inner_integral * term_factor # Note that this path.τ multiplication refer to dτ'
end


# Mexican Hat -r^2+r^4 in N-dimensions
function one_body_potential(potential::MexicanHatPotential, path::Path, bead::Int, particle::Int)
    r=norm(path.beads[mod1(bead, path.n_beads), particle,:])^2
    return 0.5 * potential.ω^2 * (-r^2+r^4)
end

# Just return value of potential for a constant potential independent of two particles.
function two_body_potential(potential::ConstantPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
    return potential.V
end

# Return the Coulomb potential between two particles.
function two_body_potential(potential::CoulombPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
    return -potential.κ / norm(path.beads[mod1(bead, path.n_beads), particle_one, :] .- path.beads[mod1(bead, path.n_beads), particle_two, :])
end

