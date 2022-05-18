# PIMC.jl

# Very generic path
mutable struct Path
	n_beads :: Int64
	n_particles :: Int64
	n_dimensions :: Int64

	beads :: Array{Float64, 3}

	τ :: Float64
	λ :: Float64

	function Path(n_beads::Int64, n_particles::Int64; n_dimensions::Int64 = 3, τ = 0.05, λ = 0.5)
		beads = rand(n_beads, n_particles, n_dimensions)
		new(n_beads, n_particles, n_dimensions, beads, τ, λ)
	end
end

function relabel_beads!(path::Path)
	rand_slice = rand(1:path.n_beads)
	slices = vcat(rand_slice:path.n_beads, 1:rand_slice-1)
	path.beads = path.beads[slices, :, :]
end

function kinetic_action(path::Path, bead::Int, particle::Int)
	bead_before = mod1(bead - 1, path.n_beads)
	bead_after = mod1(bead + 1, path.n_beads)
	kinetic_action = (norm(path.beads[bead, particle, :] - path.beads[bead_before, particle, :])^2 + norm(path.beads[bead, particle, :] - path.beads[bead_after, particle, :])^2)
	kinetic_action /= (4 * path.λ * path.τ)
	kinetic_action += 3.0 * path.n_particles / 2.0 * log(4π * path.λ * path.τ)
	return kinetic_action
end

function potential_action(path::Path, bead::Int, particle::Int, potential::ConstantPotential)
    return potential.V
end

function potential_action(path::Path, bead::Int, particle::Int, potential::OneBodyPotential)
    return path.τ * one_body_potential(potential, path, bead, particle)
end

function potential_action(path::Path, bead::Int, particle::Int, potential::TwoBodyPotential)
    potential_action = 0.0
    for particle_other in 1:path.n_particles
        if particle_other != particle
            potential_action += two_body_potential(potential, path, bead, particle, particle_other)
        end
    end
    return path.τ * potential_action
end

function potential_action(path::Path, bead::Int, particle::Int, potentials::Array{Potential})
    total_potential_action = 0.0
    for potential in potentials
        total_potential_action += potential_action(path, bead, particle, potential)
    end
    return total_potential_action
end

function primitive_action(path::Path, bead::Int, particle::Int, potential::Potential)
    return kinetic_action(path, bead, particle) + potential_action(path, bead, particle, potential)
end

function primitive_action(path::Path, bead::Int, particle::Int, potentials::Array{Potential})
    primitive_action = kinetic_action(path, bead, particle)
    for potential in potentials
        primitive_action += potential_action(path, bead, particle, potential)
    end
    return primitive_action
end

function PIMC(n_steps::Int, movers, path::Path, potentials::Union{Potential, Array{Potential}})
	observable_skip = 0.001 * n_steps
	equilibrium_skip = 0.2 * n_steps
	n_accepted = Dict(string(Symbol(mover)) => 0 for mover in movers)
	energy_trace = []
	path_trace = []
	for step in 1:n_steps

		for mover in movers
			for particle in rand(1:path.n_particles, path.n_particles)
				n_accepted[string(Symbol(mover))] += mover(path, particle, potentials)
			end
		end

		if mod(step, observable_skip) == 0 && step > equilibrium_skip
			append!(energy_trace, total_energy(path, potentials))
			append!(path_trace, [path.beads])
		end
	end

	acceptance_ratio = Dict(string(Symbol(mover)) => 1.0 * n_accepted[string(Symbol(mover))] / (n_steps * path.n_particles) for mover in movers)

	return acceptance_ratio, energy_trace, path_trace
end
