

begin
	using LinearAlgebra
	using Statistics
	using Plots
end

begin #Path construction
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
end

begin #Potential type definitions
    begin
        abstract type Potential end
        abstract type OneBodyPotential <: Potential end
        abstract type TwoBodyPotential <: Potential end
        struct ZeroPotential <: Potential end
        
        struct HarmonicPotential <: OneBodyPotential
            ω :: Float64
            function HarmonicPotential(ω::Float64)
                new(ω)
            end
        end
        
        struct FrohlichPotential <: OneBodyPotential
            α :: Float64
            ω :: Float64
            function FrohlichPotential(α::Float64, ω::Float64)
                new(α, ω)
            end
        end
        
        struct CoulombPotential <: TwoBodyPotential
            κ :: Float64
            function CoulombPotential(κ::Float64)
                new(κ)
            end
        end
    end

end

#Potential function definitons
begin
	function one_body_potential(potential::ZeroPotential, path::Path, bead::Int, particle::Int)
		return 0.0
	end

	function one_body_potential(potential::HarmonicPotential, path::Path, bead::Int, particle::Int)
		return 0.5 * potential.ω^2 * norm(path.beads[bead, particle, :])^2
	end

	function one_body_potential(potential::FrohlichPotential, path::Path, bead::Int, particle::Int)
		ω = potential.ω
		phonon_response(other_bead) = cosh(ω * path.τ * (abs(bead - other_bead) - path.n_beads / 2.0)) / sinh(ω * path.τ * path.n_beads / 2.0)
		double_integral = 0.0
		for other_bead in 1:path.n_beads
			if other_bead != bead
				double_integral += phonon_response(other_bead) / norm(path.beads[bead, particle, :] .- path.beads[other_bead, particle, :])
			end
		end
		return potential.α * ω^(3/2) / sqrt(8) * double_integral / path.n_beads
	end

	function two_body_potential(potential::ZeroPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
		return 0.0
	end

	function two_body_potential(potential::CoulombPotential, path::Path, bead::Int, particle_one::Int, particle_two::Int)
		return -potential.κ / norm(path.beads[bead, particle_one, :] .- path.beads[bead, particle_two, :])
	end
end


begin #action

    #Kinetic action
    function kinetic_action(path::Path, bead::Int, particle::Int)
        bead_before = mod1(bead - 1, path.n_beads)
        bead_after = mod1(bead + 1, path.n_beads)
        kinetic_action = (norm(path.beads[bead, particle, :] - path.beads[bead_before, particle, :])^2 + norm(path.beads[bead, particle, :] - path.beads[bead_after, particle, :])^2)
        kinetic_action /= (4 * path.λ * path.τ)
        kinetic_action += 3.0 * path.n_particles / 2.0 * log(4π * path.λ * path.τ)
        return kinetic_action
    end

    #Potential action
	function potential_action(path::Path, bead::Int, particle::Int, potential::ZeroPotential)
		return 0.0
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
end

begin #Combined action
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
end

# Energy defintion (Thermodynamic Estimator)
begin
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
end

begin #PIMC
    function PIMC(n_steps::Int, movers, observables, path::Path, potentials::Union{Potential, Array{Potential}})
        
        observable_skip = 0.001 * n_steps
        equilibrium_skip = 0.2 * n_steps
        # equilibrium_skip = 0
        
        n_accepted = Dict(string(Symbol(mover)) => 0 for mover in movers)
        observable_traces = Dict(string(Symbol(observable)) => [] for observable in observables)
        
        path_trace = []
        for step in 1:n_steps

            for mover in movers
                for particle in rand(1:path.n_particles, path.n_particles)
                    n_accepted[string(Symbol(mover))] += mover(path, particle, potentials)
                end
            end

            if mod(step, observable_skip) == 0 && step > equilibrium_skip
                for observable in observables
                    append!(observable_traces[string(Symbol(observable))], [observable(path, potentials)])
                end
                append!(path_trace, [path.beads])
            end
            
        end

        acceptance_ratio = Dict(string(Symbol(mover)) => 1.0 * n_accepted[string(Symbol(mover))] / (n_steps * path.n_particles) for mover in movers)

        return acceptance_ratio, path_trace, observable_traces
    end
end

# BELOW ARE THE DIFFERENT OBSERVABLES TO SAMPLE

begin
	function Energy(path::Path, potential::Potential)
		return kinetic_energy(path) + potential_energy(path, potential)
	end

	function Energy(path::Path, potentials::Array{Potential})
		total_energy = kinetic_energy(path)
		for potential in potentials
			total_energy += potential_energy(path, potential)
		end
		return total_energy
	end
end

begin
	function Correlation(path::Path, potential::Potential)
		correlation = Vector{Float64}(undef, path.n_beads)
		for Δτ in 1:path.n_beads, bead_one in 1:path.n_beads
			bead_two = mod1(bead_one + Δτ, path.n_beads)
			correlation[Δτ] += dot(path.beads[bead_one, :, :], path.beads[bead_two, :, :])
		end		
		return correlation ./ path.n_beads
	end
end

begin
	function Wavefunction(path::Path, potential::Potential)

	end
end

# BELOW ARE THE DIFFERENT MOVE METHODS TO SAMPLE THE BEADS

begin
    function Single(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})
        relabel_beads!(path)
        
        width = sqrt(2 * path.λ * path.τ)
        shift = width .* randn(path.n_dimensions)
        
        old_action = 0.0
        old_action += primitive_action(path, 1, particle, potentials)
        
        path.beads[1, particle, :] += shift
        
        new_action = 0.0
        new_action += primitive_action(path, 1, particle, potentials)

        if new_action - old_action <= 0.0
            return true
        elseif rand() <= exp(-(new_action - old_action))
            return true
        else
            path.beads[1, particle, :] -= shift
            return false
        end
    end


    function Displace(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})
        
        width = sqrt(2 * path.λ * path.τ)
        shift = width .* randn(path.n_dimensions)

        old_beads = copy(path.beads[:, particle, :])

        old_action = 0.0
        old_action += potential_action(path, 1, particle, potentials)

        for bead in 1:path.n_beads
            path.beads[bead, particle, :] += shift
        end

        new_action = 0.0
        new_action = potential_action(path, 1, particle, potentials)

        if new_action - old_action <= 0.0
            return true
        elseif rand() <= exp(-(new_action - old_action))
            return true
        else
            path.beads[:, particle, :] = old_beads
            return false
        end
    end


    function Staging(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})
        relabel_beads!(path)

        # segment_length = 0.025 * path.n_beads > 1 ? Int(floor(0.025 * path.n_beads)) : 1
        segment_length = 16
        
        if segment_length == 1
            return Single(path, particle, potentials)
        end

        old_beads = copy(path.beads[:, particle, :])
        old_action = 0.0
        for bead in 2:segment_length
            old_action += potential_action(path, bead, particle, potentials)
        end

        new_action = 0.0
        for bead in 1:segment_length-1
            staging_mass = (segment_length - bead + 1) / (segment_length - bead)
            staging_position = (path.beads[1 + segment_length, particle, :] + path.beads[bead, particle, :] * (segment_length - bead)) / (segment_length - bead + 1)
            path.beads[bead + 1, particle, :] = staging_position + randn(path.n_dimensions) * sqrt(path.τ / staging_mass)
            new_action += potential_action(path, bead + 1, particle, potentials)
        end

        if new_action - old_action < 0.0
            return true
        elseif rand() < exp(-(new_action - old_action))
            return true
        else
            path.beads[:, particle, :] = old_beads
            return false
        end
    end

    function Bisection(path::Path, particle::Int, potentials::Union{Potential, Array{Potential}})
        relabel_beads!(path)

        max_level = Int(floor(log(rand(1:path.n_beads)) / log(2)))
        # max_level = 3
        clip_length = 2^max_level + 1

        old_beads = copy(path.beads[:, particle, :])
        
        total_old_action = 0.0
        for bead in 2:clip_length-1
            total_old_action += potential_action(path, bead, particle, potentials)
        end

        old_action = 0.0
        new_action = 0.0
        for level in max_level:-1:1
            step = 2^(level - 1)
            ratio = 2^max_level / step
            for interval in 1:2:ratio
                bead = Int(1 + interval * 2^max_level / ratio)
                old_action += potential_action(path, bead, particle, potentials)
                shift_vector = randn(path.n_dimensions) * sqrt(step * path.τ * path.λ)
                path.beads[bead, particle, :] = 0.5 * (path.beads[bead - step, particle, :] + path.beads[bead + step, particle, :]) + shift_vector
                new_action += potential_action(path, bead, particle, potentials)
            end
            if rand() >= exp(-(new_action - old_action))
                return false
            end
        end
        
        total_new_action = 0.0
        for bead in 2:clip_length-1
            total_new_action += potential_action(path, bead, particle, potentials)
        end

        if total_new_action - total_old_action < 0.0
            return true
        elseif rand() < exp(-(total_new_action - total_old_action))
            return true
        else
            path.beads[:, particle, :] = old_beads
            return false
        end
    end

end

# ╔═╡ 8d8a2113-a6f7-4c40-b04e-8e8ca16a70b3
function draw_beads_3d(path, xlims, ylims, zlims)

	p = plot()
	
	for particle in 1:size(path)[2]
		x = path[:, particle, 1]
		y = path[:, particle, 2]
		z = path[:, particle, 3]
		
		x = reshape(x, length(x))
		y = reshape(y, length(y))
		z = reshape(z, length(z))

		push!(x, x[1])
		push!(y, y[1])
		push!(z, z[1])

		plot!(p, x, y, z, marker = :circle, label = "P $particle", legend = false, xlims = xlims, ylims = ylims, zlims = zlims)
	end
	return p
end

#actual PIMC Code
begin
	T = 1.0
	λ = 0.5
	n_beads = 100
	τ = 1.0 / (T * n_beads)
	n_particles = 1
	path = Path(n_beads, n_particles, τ = τ, λ = λ)

	n_steps = 1000000
	p_beads = copy(path.beads)
	pimc = PIMC(n_steps, [Single, Displace], [Energy, Correlation], path, HarmonicPotential(1.0))
end

# ╔═╡ 2c80c5f3-f09e-4c95-9667-f10b99260f92
equipartition_energy = 3/2 * T * n_particles

# ╔═╡ 396be977-6266-4f1e-a9f9-bd1d3c03e6b7
harmonic_potential_energy = 3.0 / 2.0 * coth(1.0 / 2.0 / T)

