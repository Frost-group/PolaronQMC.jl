# actions.jl

function kinetic_action(path::Path, bead_one::Int, bead_two::Int, particle::Int)
	kinetic_action = norm(path.beads[bead_two, particle, :] - path.beads[bead_one, particle, :])^2
	kinetic_action /= (4 * path.λ * path.τ)
	kinetic_action += 3.0 * path.n_particles * path.n_beads / 2.0 * log(4π * path.λ * path.τ)
	return kinetic_action
end

function potential_action(path::Path, bead::Int, particle::Int, potential::ConstantPotential)
    return path.τ * potential.V
end

function potential_action(path::Path, bead::Int, particle::Int, potential::OneBodyPotential)
    return path.τ * one_body_potential(potential, path, bead, particle)
end

function potential_action(path::Path, bead::Int, particle::Int, potential::TwoBodyPotential)
    potential_action = sum(two_body_potential(potential, path, bead, particle, other_particle) for other_particle in 1:path.n_particles if particle != other_particle)
    return path.τ * potential_action
end

function potential_action(path::Path, bead::Int, particle::Int, potentials::Array{Potential})
    sum(potential_action(path, bead, particle, potential) for potential in potentials)
end

function primitive_action(path::Path, bead_one::Int, bead_two::Int, particle::Int, potential::Potential)
    return kinetic_action(path, bead_one, bead_two, particle) + potential_action(path, bead_two, particle, potential)
end

function primitive_action(path::Path, bead_one::Int, bead_two::Int, particle::Int, potentials::Array{Potential})
    sum(kinetic_action(path, bead_one, bead_two, particle), potential_action(path, bead_two, particle, potential) for potential in potentials)
end
