# actions.jl
# New comment testing
# New commit 2

#Ways of calculating action

    #simple method


    function kinetic_action(path::Path, bead_one::Int, bead_two::Int, particle::Int, regime::Simple_Regime)
        kinetic_action = 0.5 * path.m * norm(path.beads[bead_two, particle, :] - path.beads[bead_one, particle, :])^2
        return kinetic_action
    end


    function potential_action(path::Path, bead::Int, particle::Int, potential::OneBodyPotential, regime::Simple_Regime)
        return one_body_potential(potential, path, bead, particle)
    end


    function total_action(path::Path, bead_one::Int, bead_two::Int, particle::Int, potential::OneBodyPotential, regime::Simple_Regime)
        return kinectic_action(path, bead_one, bead_two, particle, regime) + potential_action(path, bead_one, particle, potential, regime)
    end



#Primitive method (based off Ceperly paper)

"""
    kinetic_action(path::Path, bead_one::Int, bead_two::Int, particle::Int)

Kinetic action for a link between two beads `bead_one` and `bead_two` for a particle indexed by an integer `particle`. Usually, these are neighbouring beads.

# Arguments
- `path::Path`: collection of all particle imaginary-time paths.
- `bead_one::Int`: first bead in the link.
- `bead_two::Int`: second bead in the link.
- `particle::Int`: select a specific particle indexed by this integer.

See also [`Path`](@ref). 
"""
function kinetic_action(path::Path, bead_one::Int, bead_two::Int, particle::Int, regime::Primitive_Regime)

    # Contribution from per degree of freedom per particle for distinguishable particles.
    # Comes from the normalisation term of the density matrix. Can be ignored but is included for completeness.
    kinetic_action = path.n_dimensions * path.n_particles / 2.0 * log(4π * path.λ * path.τ)

    # Contribution from the link connecting the two beads.
	kinetic_action += norm(path.beads[bead_two, particle, :] - path.beads[bead_one, particle, :])^2 / (4 * path.λ * path.τ)

	return kinetic_action
end


"""
    potential_action(path::Path, bead::Int, particle::Int, potential::ConstantPotential)

Potential action for a bead indexed by an integer `bead` and a particle indexed by an integer `particle`.
This is the potential action for an external constant potential.

# Arguments
- `path::Path`: collection of all particle imaginary-time paths.
- `bead::Int`: select a specific bead indexed by an integer.
- `particle::Int`: select a specific particle indexed by this integer.
- `potentials::ConstantPotential`: a constant potential type. Typically just zero.

See also [`Path`](@ref), [`ConstantPotential`](@ref). 
"""
function potential_action(path::Path, bead::Int, particle::Int, potential::ConstantPotential, regime::Primitive_Regime)
    return path.τ * potential.V
end


"""
    potential_action(path::Path, bead::Int, particle::Int, potential::OneBodyPotential)

Potential action for a bead indexed by an integer `bead` and a particle indexed by an integer `particle`.
This is the potential action for an external one-body potential.

# Arguments
- `path::Path`: collection of all particle imaginary-time paths.
- `bead::Int`: select a specific bead indexed by an integer.
- `particle::Int`: select a specific particle indexed by this integer.
- `potentials::OneBodyPotential`: a one-body potential type. Acts like an external potential on the system.

See also [`Path`](@ref), [`OneBodyPotential`](@ref). 
"""
function potential_action(path::Path, bead::Int, particle::Int, potential::OneBodyPotential, regime::Primitive_Regime)
    return path.τ * one_body_potential(potential, path, bead, particle)
end


"""
    potential_action(path::Path, bead::Int, particle::Int, potential::TwoBodyPotential)

Potential action for a bead indexed by an integer `bead` and a particle indexed by an integer `particle`.
This is the potential action for an internal two-body potential between the specified bead on the specified particle and the same bead on all other particles. 

# Arguments
- `path::Path`: collection of all particle imaginary-time paths.
- `bead::Int`: select a specific bead indexed by an integer.
- `particle::Int`: select a specific particle indexed by this integer.
- `potentials::TwoBodyPotential`: a two-body potential type. Acts like interactions between pairs of particles in the system.

See also [`Path`](@ref), [`TwoBodyPotential`](@ref). 
"""
function potential_action(path::Path, bead::Int, particle::Int, potential::TwoBodyPotential, regime::Primitive_Regime)
    potential_action = sum(two_body_potential(potential, path, bead, particle, other_particle) for other_particle in 1:path.n_particles if particle != other_particle)
    return path.τ * potential_action 
end


function total_action(path::Path, bead_one::Int, bead_two::Int, particle::Int, potential::Potential, regime::Primitive_Regime)
        return kinetic_action(path, bead_one, bead_two, particle, regime) + potential_action(path, bead_two, particle, potential, regime)
end
