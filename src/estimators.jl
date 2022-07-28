# estimators.jl

function kinetic_energy(path::Path, potentials)
    bead = rand(1:path.n_beads)
	prefactor = -1.0 / (4.0 * path.λ * path.τ^2)
	kinetic_energy = sum(norm(path.beads[bead, particle, :] - path.beads[bead-1, particle, :])^2 for particle in 1:path.n_particles)
	return prefactor * kinetic_energy + path.n_dimensions * path.n_particles / (2 * path.τ)
end

function potential_energy(path::Path, potential::ConstantPotential)
    return potential.V
end

function potential_energy(path::Path, potential::OneBodyPotential)
    bead = rand(1:path.n_beads)
    potential_energy = sum(one_body_potential(potential, path, bead, particle) for particle in 1:path.n_particles)
    return potential_energy
end

function potential_energy(path::Path, potential::TwoBodyPotential)
    bead = rand(1:path.n_beads)
    potential_energy = sum(
        two_body_potential(potential, path, bead, particle_one, particle_two) 
        for particle_one in 1:path.n_particles, 
            particle_two in 1:(particle_one-1)
            )
    return potential_energy
end

function potential_energy(path::Path, potentials::Array{Potential})
    potential_energy = sum(potential_energy.(path, potentials))
    return potential_energy
end

function Energy(path::Path, potential::Potential)
    return kinetic_energy(path, potential) + potential_energy(path, potential)
end

function Energy(path::Path, potentials::Array{Potential})
    total_energy = kinetic_energy(path, potentials)
    total_energy += sum(potential_energy.(path, potentials))
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

