# observables.jl

function kinetic_energy(path::Path)
	kinetic_energy = 0.0
	prefactor = 1.0 / (4.0 * path.λ * path.τ^2)
	for bead_one in 1:path.n_beads
		bead_two = mod1(bead_one - 1, path.n_beads)
		for particle in 1:path.n_particles
			kinetic_energy -= norm(path.beads[bead_one, particle, :] - path.beads[bead_two, particle, :])^2
		end
	end
	return prefactor * kinetic_energy / path.n_beads + path.n_dimensions * path.n_particles / (2 * path.τ)
end

function potential_energy(path::Path, potential::ZeroPotential)
    return 0.0
end

function potential_energy(path::Path, potential::OneBodyPotential)
    potential_energy = 0.0
    for bead in 1:path.n_beads, particle in 1:path.n_particles
        potential_energy += one_body_potential(potential, path, bead, particle)
    end
    return potential_energy / path.n_beads
end

function potential_energy(path::Path, potential::TwoBodyPotential)
    potential_energy = 0.0
    for bead in 1:path.n_beads, particle_one in 1:path.n_particles, particle_two in 1:path.n_particles
        if particle_one != particle_two
            potential_energy += two_body_potential(potential, path, bead, particle_one, particle_two)
        end
    end
    return potential_energy / path.n_beads
end

function potential_energy(path::Path, potentials::Array{Potential})
    total_potential_energy = 0.0
    for potential in potentials
        total_potential_energy += potential_energy(path, potential)
    end
    return total_potential_energy
end
