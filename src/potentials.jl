# potentials.jl

"""
Collection of potentials to use in simulations. 
"""

"""
Outer constructors for different potential types.
"""
# Defining delta function

# Just return value of potential for a constant potential independent of single particle.
function oneBodyPotential(
    potential::ConstantPotential,
    path::Path,
    bead::Int,
    particle::Int,
)
    return potential.V
end

# Return the harmonic potential for a single particle.
function oneBodyPotential(
    potential::HarmonicPotential,
    path::Path,
    bead::Int,
    particle::Int,
    store_diff::Vector{Float64},
    prop_Matrix::Array{Float64},
)
    return 0.5 *
           path.m *
           potential.ω^2 *
           norm(path.beads[mod1(bead, path.n_beads), particle, :])^2
end

function oneBodyPotential(
    potential::HarmonicInteractionPotential,
    path::Path,
    bead::Int,
    particle::Int,
)
    harmonic_pot =
        0.5 *
        path.m *
        potential.ω^2 *
        norm(path.beads[mod1(bead, path.n_beads), particle, :])^2
    coloumb_int = 0
    for other_particle = 1:path.n_particles
        if other_particle != particle
            coloumb_int +=
                -potential.κ / norm(
                    path.beads[mod1(bead, path.n_beads), particle, :] .-
                    path.beads[mod1(bead, path.n_beads), other_particle, :],
                )
        end
    end
    return harmonic_pot + coloumb_int
end

# Returns the Frohlich potential for a single particle
function oneBodyPotential(
    potential::FrohlichPotential,
    path::Path,
    bead::Int,
    particle::Int,
    store_diff::Vector{Float64},
    prop_Matrix::Array{Float64},
)
    # Refer to the Lecture notes on Frohlich Polaron by Devreese
    # Allowing sum over multi-phonon modes

    # Defining the constants to avoid repeated attributes call
    n_beads, l = path.n_beads, length(potential.ω)
    β, α = path.τ * n_beads, potential.α
    #term_factor = -0.5 * α * (ħω)^(3/2) * sqrt(1/2/path.m) * csch(ħω * β / 2);

    # Calculates the double integral component with individual modes
    final_integral = 0.0 # The final after summing all beads
    for other_bead = 1:n_beads
        if other_bead != bead # Discounting self-contributions
            for i = 1:path.n_dimensions
                store_diff[i] =
                    path.beads[bead, particle, i] -
                    path.beads[mod1(other_bead, n_beads), particle, i]
            end

            for i = 1:l
                final_integral +=
                    0.5 *
                    α *
                    (potential.ħ * potential.ω[i])^(3 / 2) *
                    sqrt(1 / 2 / path.m) *
                    csch(potential.ħ * potential.ω[i] * β / 2) *
                    prop_Matrix[
                        max(bead, mod1(other_bead, n_beads)),
                        min(bead, mod1(other_bead, n_beads)),
                        i,
                    ] / norm(store_diff)
            end
            #A = 1/(norm(path.beads[mod1(bead, path.n_beads), particle, :] - path.beads[mod1(other_bead, path.n_beads), particle, :]))
            #inner_integral += prop_Matrix[max(bead, other_bead), min(bead, other_bead)]/norm(store_diff)
        end
    end
    return -path.τ * final_integral #* term_factor # Note that this path.τ multiplication refer to dτ', the minus sign is due to negative potential
end

function oneBodyPotential(
    potential::FrohlichInteractionPotential,
    path::Path,
    bead::Int,
    particle::Int,
)
    β = path.τ * path.n_beads
    m = path.m
    ħω = potential.ħ * potential.ω
    α = potential.α
    term_factor = -0.5 * α * (ħω)^(3 / 2) * sqrt(1 / 2 / m) * csch(ħω * β / 2)

    inner_integral = 0.0
    for other_bead = 1:path.n_beads
        if other_bead != bead
            inner_integral +=
                cosh(ħω * β * (abs(bead - other_bead) / path.n_beads - 0.5)) / norm(
                    path.beads[mod1(bead, path.n_beads), particle, :] -
                    path.beads[mod1(other_bead, path.n_beads), particle, :],
                )
        end
    end

    coloumb_int = 0
    for other_particle = 1:path.n_particles
        if other_particle != particle
            coloumb_int +=
                -potential.κ / norm(
                    path.beads[mod1(bead, path.n_beads), particle, :] .-
                    path.beads[mod1(bead, path.n_beads), other_particle, :],
                )
        end
    end
    return path.τ * inner_integral * term_factor + coloumb_int # Note that this path.τ multiplication refer to dτ'
end

# Mexican Hat -r^2+r^4 in N-dimensions
function oneBodyPotential(
    potential::MexicanHatPotential,
    path::Path,
    bead::Int,
    particle::Int,
)
    r = norm(path.beads[mod1(bead, path.n_beads), particle, :])^2
    return 0.5 * potential.ω^2 * (-r^2 + r^4)
end

function oneBodyPotential(
    potential::HolsteinPotential,
    path::DiscretePath,
    bead::Int,
    particle::Int64,
    F_l::Array{Float64},
)
    # Refer to the Holstein Small-polaron paper

    # Defining the constants to avoid repeated attributes call
    final_integral = 0.0
    for other_bead = 1:path.n_beads
        if bead != other_bead
            if δ(
                @view(path.beads[bead, particle, :]),
                @view(path.beads[other_bead, particle, :])
            )
                #if @views path.beads[bead, particle, :] != @views path.beads[other_bead, particle, :]
                final_integral += F_l[abs(bead - other_bead)+1]
            end
        end
    end
    return final_integral #* term_factor # Note that this path.τ multiplication refer to dτ'
end


# Just return value of potential for a constant potential independent of two particles.
function twoBodyPotential(
    potential::ConstantPotential,
    path::Path,
    bead::Int,
    particle_one::Int,
    particle_two::Int,
)
    return potential.V
end


# Return the Coulomb potential between two particles.
function twoBodyPotential(
    potential::CoulombPotential,
    path::Path,
    bead::Int,
    particle_one::Int,
    particle_two::Int,
)
    return -potential.κ / norm(
        path.beads[mod1(bead, path.n_beads), particle_one, :] .-
        path.beads[mod1(bead, path.n_beads), particle_two, :],
    )
end
