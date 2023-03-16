# actions.jl


# Ways of calculating action

# not used?
function totalAction(path::Path, bead_one::Int, bead_two::Int, particle::Int, potential::OneBodyPotential, regime::Regime)
    return kinecticAction(path, bead_one, bead_two, particle, regime) + potential_action(path, bead_one, particle, potential, regime)
end

"""
Calculate Kinetic Actions
"""


function kineticAction(path::Path, bead_one::Int, bead_two::Int, particle::Int, regime::SimpleRegime)
    kinetic_action = 0.5 * path.m * norm(path.beads[mod1(bead_two, path.n_beads), particle, :] - path.beads[mod1(bead_one, path.n_beads), particle, :])^2 / path.τ
    return kinetic_action
end

function kineticAction(path::Path, bead_one::Int, bead_two::Int, particle::Int, regime::Regime)
    kinetic_action = 0.5 * path.m * norm(path.beads[mod1(bead_two, path.n_beads), particle, :] - path.beads[mod1(bead_one, path.n_beads), particle, :])^2 / path.τ
    return kinetic_action
end

#Primitive method (based off Ceperly paper)
function kineticAction(path::Path, bead_one::Int, bead_two::Int, particle::Int, regime::PrimitiveRegime)

    """
    kinetic_action(path::Path, bead_one::Int, bead_two::Int, particle::Int)

    Kinetic action for a link between two beads `bead_one` and `bead_two` for a particle indexed by an integer `particle`. Usually, these are neighbouring beads.

    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead_one::Int`: first bead in the link.
    - `bead_two::Int`: second bead in the link.
    - `particle::Int`: select a specific particle indexed by this integer.

    See also [`Path`](@ref). 
    """

    # Contribution from per degree of freedom per particle for distinguishable particles.
    # Comes from the normalisation term of the density matrix. Can be ignored but is included for completeness.
    #kinetic_action = path.n_dimensions * path.n_particles / 2.0 * log(4π * path.λ * path.τ)

    # Contribution from the link connecting the two beads.
	#kinetic_action += norm(path.beads[mod1(bead_two, path.n_beads), particle, :] - path.beads[mod1(bead_one, path.n_beads), particle, :])^2 / (4 * path.λ * path.τ)
    kinetic_action = norm(path.beads[mod1(bead_two, path.n_beads), particle, :] - path.beads[mod1(bead_one, path.n_beads), particle, :])^2 / (4 * path.λ * path.τ)
	return kinetic_action
end

"""
Calculate Potential Action
"""

function potentialAction(path::Path, bead::Int, particle::Int, potential::OneBodyPotential, regime::SimpleRegime)
    return oneBodyPotential(potential, path, bead, particle) * path.τ
end


function potentialAction(path::Path, bead::Int, particle::Int, potential::ConstantPotential, regime::PrimitiveRegime)

    """
    potential_action(path::Path, bead::Int, particle::Int, potential::ConstantPotential)

    Potential action for a bead indexed by an integer `bead` and a particle indexed by an integer `particle`.
    This is the potential action for an external constant potential.

    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead::Int`: select a specific bead indexed by an integer.
    - `particle::Int`: select a specific particle indexed by this integer.
    - `potentials::ConstantPotential`: a constant potential type. Typically just zero.

    See also [`Path`](@ref), [`ConstantPotential`](@ref). 
    """

    return path.τ * potential.V
end

function potentialAction(path::Path, bead::Int, particle::Int, potential::FrohlichPotential, regime::LBRegime)
    β = path.τ * path.n_beads
    t2_prefactor = path.τ^3 / (24 * path.m)
    V2_prefactor = (0.5 * potential.α * (potential.ħ * potential.ω)^(3/2) * sqrt(1/2/path.m) * csch(potential.ħ * potential.ω * β / 2))^2
    ħω = potential.ω * potential.ħ
    # F = -dV/dr ∝ r/(r(τ)-r(τ'))^3
    
    inner_integral = 0.0
    
    #if bead != other_bead
    for other_bead in 1:path.n_beads+1
        if mod1(bead, path.n_beads) != mod1(other_bead, path.n_beads)
            g_factor = (cosh(ħω * β * (abs(bead-other_bead)/(path.n_beads) - 0.5)))^2
            inner_integral += g_factor / norm(path.beads[mod1(bead, path.n_beads),particle,:] - path.beads[mod1(other_bead, path.n_beads),particle,:])^4
        end
    end
    return 2 * path.τ * oneBodyPotential(potential, path, bead, particle) + 2 * inner_integral * path.τ^2 * t2_prefactor * V2_prefactor
end


function potentialAction(path::Path, bead::Int, particle::Int, potential::FrohlichPotential, regime::BoundRegime)
    n_beads = path.n_beads
    β = path.τ * path.n_beads
    m = path.m
    ħω = potential.ħ * potential.ω
    α = potential.α
    term_factor = -0.5 * α * (ħω)^(3/2) * sqrt(1/2/m) * csch(ħω * β / 2)
    
    inner_integral = 0.0
    special_integral = 0.0
    if mod1(bead, n_beads) == 1
        for other_bead in 2:n_beads
            special_integral += cosh(ħω * β * (abs(1-other_bead)/n_beads - 0.5)) / norm(path.beads[mod1(bead, n_beads), particle, :] - path.beads[mod1(other_bead, n_beads), particle, :])
            special_integral += cosh(ħω * β * (abs(n_beads+1-other_bead)/n_beads - 0.5)) / norm(path.beads[mod1(bead, n_beads), particle, :] - path.beads[mod1(other_bead, n_beads), particle, :])
        end
    else
        for other_bead in 1:n_beads+1
            if mod1(other_bead, n_beads) == 1
                special_integral += cosh(ħω * β * (abs(bead-other_bead)/n_beads - 0.5)) / norm(path.beads[mod1(bead, n_beads), particle, :] - path.beads[mod1(other_bead, n_beads), particle, :])
            else
                if mod1(other_bead, n_beads) != mod1(bead, n_beads)
                    inner_integral += cosh(ħω * β * (abs(bead-other_bead)/n_beads - 0.5)) / norm(path.beads[mod1(bead, n_beads), particle, :] - path.beads[mod1(other_bead, n_beads), particle, :])
                    #inner_integral += g_factor / norm(path.beads[mod1(bead, path.n_beads), particle, :] - path.beads[mod1(other_bead, path.n_beads), particle, :])
                end
            end
        end
    end
    return (2 * inner_integral + special_integral) * term_factor * path.τ^2 # Note that this path.τ multiplication refer to dτ'
end

function bisectPotentialAction(path::Path, bead::Int, beadrange::Union{Vector{Int64}, StepRange{Int64, Int64}, UnitRange{Int64}}, particle::Int, potential::FrohlichPotential, regime::PrimitiveRegime)
    β = path.τ * path.n_beads
    m = path.m
    ħω = potential.ħ * potential.ω
    α = potential.α
    term_factor = -0.5 * α * (ħω)^(3/2) * sqrt(1/2/m) * csch(ħω * β / 2)
    n_beads = path.n_beads
    bead = mod1(bead, n_beads)
    inner_integral = 0.0

    for other_bead in 1:path.n_beads
        if other_bead != bead
            if other_bead in beadrange
                inner_integral += 0.5 * cosh(ħω * β * (abs(bead-other_bead)/n_beads - 0.5)) / norm(path.beads[mod1(bead, n_beads), particle, :] - path.beads[mod1(other_bead, n_beads), particle, :])
            else
                inner_integral += cosh(ħω * β * (abs(bead-other_bead)/n_beads - 0.5)) / norm(path.beads[mod1(bead, n_beads), particle, :] - path.beads[mod1(other_bead, n_beads), particle, :])
            end
        end
    end
    return 2 * (path.τ)^2 * inner_integral * term_factor # Note that this path.τ multiplication refer to dτ'
end

function potentialAction(path::Path, bead::Int, particle::Int, potential::OneBodyPotential, regime::PrimitiveRegime)

    """
    potential_action(path::Path, bead::Int, particle::Int, potential::OneBodyPotential)

    Potential action for a bead indexed by an integer `bead` and a particle indexed by an integer `particle`.
    This is the potential action for an external one-body potential.

    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead::Int`: select a specific bead indexed by an integer.
    - `particle::Int`: select a specific particle indexed by this integer.
    - `potentials::OneBodyPotential`: a one-body potential type. Acts like an external potential on the system.

    See also [`Path`](@ref), [`OneBodyPotential`](@ref). 
    """
    
    if typeof(potential) == FrohlichPotential
        return 2 * path.τ * oneBodyPotential(potential, path, bead, particle)
    end

    return path.τ * oneBodyPotential(potential, path, bead, particle)
end


function potentialAction(path::Path, bead::Int, particle::Int, potential::TwoBodyPotential, regime::PrimitiveRegime)

    """
    potential_action(path::Path, bead::Int, particle::Int, potential::TwoBodyPotential)

    Potential action for a bead indexed by an integer `bead` and a particle indexed by an integer `particle`.
    This is the potential action for an internal two-body potential between the specified bead on the specified particle and the same bead on all other particles. 

    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead::Int`: select a specific bead indexed by an integer.
    - `particle::Int`: select a specific particle indexed by this integer.
    - `potentials::TwoBodyPotential`: a two-body potential type. Acts like interactions between pairs of particles in the system.

    See also [`Path`](@ref), [`TwoBodyPotential`](@ref). 
    """

    potential_action = sum(twoBodyPotential(potential, path, bead, particle, other_particle) for other_particle in 1:path.n_particles if particle != other_particle)
    return path.τ * potential_action 
end
