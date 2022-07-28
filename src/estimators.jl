# estimators.jl

#types of estimators
abstract type Estimator end

#Energy estimators

struct Simple_Estimator <: Estimator end #Estimator using basic sum of kinetic and potential total_energy

struct Simple_Virial_Estimator <: Estimator end #Estimator using virial theorem for potential term

struct Thermodynamic_Estimator <: Estimator end #Estimator using thermodynamic theory

struct Virial_Estimator <: Estimator
    L :: Int
    function Virial_Estimator(L::Int)
        new(L)
    end
end




#Kinectic energy estimators ---------------------------------

function kinetic_energy(path::Path, estimator::Union{Simple_Estimator, Simple_Virial_Estimator}) #thermal dynamic estimator from ceperly paper
	kinetic_energy = 0.0
	for bead in 1:path.n_beads, particle in 1:path.n_particles
		kinetic_energy += 0.5*path.m*(path.beads[bead, particle] - path.beads[bead-1, particle])^2
	end
	return kinetic_energy / path.n_beads
end



function kinetic_energy(path::Path, estimator::Thermodynamic_Estimator) #thermal dynamic estimator from ceperly paper
	
    kinetic_energy = 0.0
    term_one = path.n_particles / (2*path.τ)
	prefactor = 1.0 / (4.0 * path.λ * path.τ^2)
	for bead in 1:path.n_beads, particle in 1:path.n_particles
		kinetic_energy += (path.beads[bead, particle] - path.beads[bead-1, particle])^2
	end
    return term_one - (kinetic_energy * prefactor/ path.n_beads)
end

function kinetic_energy(path::Path, estimator::Virial_Estimator)
    kinetic_energy = 0.0
    term_one = 1.0/(2 * estimator.L * path.τ)
    prefactor = 1.0 / (4.0 * estimator.L * path.τ^2 * path.λ)
    for bead in 1:path.n_beads, particle in 1:path.n_particles
        link_term = path.beads[bead+estimator.L, particle] - path.beads[bead, particle]
        link_term *= path.beads[bead+1, particle] - path.beads[bead, particle]
        kinetic_energy += (term_one - prefactor * link_term)
    end
    return kinetic_energy / path.n_beads
end


#Potential energy estimators -------------------------------------------

function potential_energy(path::Path, potential::OneBodyPotential, estimator::Union{Simple_Estimator, Simple_Virial_Estimator, Thermodynamic_Estimator})
    potential_energy = 0.0
    for bead in 1:path.n_beads, particle in 1:path.n_particles
        potential_energy += one_body_potential(potential, path, bead, particle)
    end
    return potential_energy / path.n_beads
end



function potential_energy(path::Path, potential::ConstantPotential, estimator::Thermodynamic_Estimator)
    return potential.V
end


function potential_energy(path::Path, potential::HarmonicPotential, estimator::Virial_Estimator)
    potential_energy = 0.0
    for bead in 1:path.n_beads, particle in 1:path.n_particles
        Δ = 0.0

        for j in -1*estimator.L+1:estimator.L-1
            Δ += path.beads[bead, particle] - path.beads[bead+j, particle]
        end
        Δ /= 2*estimator.L

        generalised_force = -1 / path.τ * -path.m * potential.ω^2 * path.beads[bead, particle]
        potential_energy += -0.5*generalised_force*Δ + one_body_potential(potential, path, bead, particle)
    end
    return potential_energy / path.n_beads
end

function potential_energy(path::Path, potential::TwoBodyPotential, estimator::Thermodynamic_Estimator)
    potential_energy = 0.0
    for bead in 1:path.n_beads, particle_one in 1:path.n_particles, particle_two in 1:path.n_particles
        if particle_one != particle_two
            potential_energy += two_body_potential(potential, path, bead, particle_one, particle_two)
        end
    end
    return potential_energy / path.n_beads
end


# Energy ---------------------------------------------------------



function Energy(path::Path, potential::Potential, estimator::Estimator)
    return kinetic_energy(path, estimator) + potential_energy(path, potential, estimator)
end


# Correlation ---------------------------------------------------------------------

function Correlation(path::Path, potential::Potential, estimator::Estimator)
    correlation = Vector{Float64}(undef, path.n_beads)
    for Δτ in 1:path.n_beads, bead_one in 1:path.n_beads
        bead_two = mod1(bead_one + Δτ, path.n_beads)
        correlation[Δτ] += dot(path.beads[bead_one, :, :], path.beads[bead_two, :, :])
    end		
    return correlation ./ path.n_beads
end


