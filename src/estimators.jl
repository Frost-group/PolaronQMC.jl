# estimators.jl

function kinetic_energy(path::Path, potentials)
	kinetic_energy = 0.0
    bead = rand(1:path.n_beads)
	prefactor = -1.0 / (4.0 * path.λ * path.τ^2)
	for particle in 1:path.n_particles
		kinetic_energy += norm(path.beads[bead, particle, :] - path.beads[bead-1, particle, :])^2
	end
	return prefactor * kinetic_energy + path.n_dimensions * path.n_particles / (2 * path.τ)
end

function potential_energy(path::Path, potential::ConstantPotential)
    return potential.V
end

function potential_energy(path::Path, potential::OneBodyPotential)
    potential_energy = 0.0
    bead = rand(1:path.n_beads)
    for particle in 1:path.n_particles
        potential_energy += one_body_potential(potential, path, bead, particle)
    end

    # potential_energy = primitive_inter_action(path, bead, bead-1, 1, potential) / path.τ
    return potential_energy
end

function potential_energy(path::Path, potential::TwoBodyPotential)
    potential_energy = 0.0
    bead = rand(1:path.n_beads)
    for particle_one in 1:path.n_particles, particle_two in 1:path.n_particles
        if particle_one != particle_two
            potential_energy += two_body_potential(potential, path, bead, particle_one, particle_two)
        end
    end
    return potential_energy
end

function potential_energy(path::Path, potentials::Array{Potential})
    total_potential_energy = 0.0
    for potential in potentials
        total_potential_energy += potential_energy(path, potential)
    end
    return total_potential_energy
end

function Energy(path::Path, potential::Potential)
    return kinetic_energy(path, potential) + potential_energy(path, potential)
end

function Energy(path::Path, potentials::Array{Potential})
    total_energy = kinetic_energy(path, potentials)
    for potential in potentials
        total_energy += potential_energy(path, potential)
    end
    return total_energy
end

# function Virial_Energy(path::Path, potential::Potential; window_size = 1)
#     window_size = path.n_beads
#     total_energy = path.n_dimensions * path.n_particles / (2 * window_size * path.τ * path.n_beads)
#     general_force = 0.0
#     for bead in 1:path.n_beads, particle in 1:path.n_particles
#         general_force -= primit
#     end

# end

function Correlation(path::Path, potential::Potential)
    correlation = Vector{Float64}(undef, path.n_beads)
    for Δτ in 1:path.n_beads, bead_one in 1:path.n_beads
        bead_two = mod1(bead_one + Δτ, path.n_beads)
        correlation[Δτ] += dot(path.beads[bead_one, :, :], path.beads[bead_two, :, :])
    end		
    return correlation ./ path.n_beads
end

