# estimators.jl
using Statistics
using ThreadsX

#types of estimators
abstract type Estimator end

#Energy estimators

struct Simple_Estimator <: Estimator end #Estimator using basic sum of kinetic and potential total_energy

struct Simple_Virial_Estimator <: Estimator end #Estimator using virial theorem for potential term

struct Thermodynamic_Estimator <: Estimator end #Estimator using thermodynamic theory

struct Virial_Estimator <: Estimator end #Estimator derived using virial theorem

struct Virial_EstimatorX <: Estimator end #Estimator derived using virial theorem




#Kinectic energy estimators ---------------------------------

function kinetic_energy(path::Path, potential::Potential, estimator::Union{Simple_Estimator, Simple_Virial_Estimator}) #thermal dynamic estimator from ceperly paper
	kinetic_energy = 0.0
	for bead in 1:path.n_beads, particle in 1:path.n_particles
		kinetic_energy += 0.5*path.m*(path.beads[bead, particle] - path.beads[bead-1, particle])^2
	end
	return kinetic_energy / path.n_beads
end



function kinetic_energy(path::Path, potential::Potential, estimator::Thermodynamic_Estimator) #thermal dynamic estimator from ceperly paper
	
    kinetic_energy = 0.0
    term_one = path.n_dimensions * path.n_particles / (2*path.τ)
	prefactor = 1.0 / (4.0 * path.λ * path.τ^2)
	for bead in 1:path.n_beads, particle in 1:path.n_particles
		kinetic_energy += norm(path.beads[bead, particle, :] - path.beads[bead-1, particle, :])^2
	end
    return term_one - (kinetic_energy * prefactor / path.n_beads)
end




function kinetic_energy(path::Path, potential::HarmonicPotential, estimator::Virial_Estimator)
    term_one = (path.n_dimensions * path.n_particles) / (2 * path.τ * path.n_beads) 

    #term 2 prefactor
    t2_prefactor = (path.m * potential.ω^2) / (2 * path.n_beads)
    term_two = 0.0
    for particle in 1:path.n_particles
        centroid_pos = zeros(path.n_dimensions)
        for bead in 1:path.n_beads
            centroid_pos += path.beads[bead,particle,:]
        end

        for bead in 1:path.n_beads
            term_two += dot(path.beads[bead, particle, :] - centroid_pos, path.beads[bead, particle, :])
        end
    end

    return term_one + t2_prefactor * term_two

end



function kinetic_energy(path::Path, potential::FrohlichPotential, estimator::Virial_Estimator)
    term_one = path.n_dimensions * path.n_particles / (2 * path.τ * path.n_beads)
    #term prefactor
    t2_prefactor = path.τ / (2 * path.n_beads)

    #beta conversion g_factor
    #β_conversion_factor = 47.81
    β_conversion_factor = 1

    β_reduced = path.τ * path.n_beads * β_conversion_factor


    function get_term_two(particle, bead, other_bead, centroid_pos)
        if bead != other_bead
            g_factor = 0.5 * potential.α * sqrt(potential.ħ/(2*path.m*potential.ω)) * cosh(β_reduced  * (abs(bead - other_bead)/path.n_beads - 0.5)) * csch(β_reduced * 0.5)
            return g_factor * dot((path.beads[bead,particle,:] - centroid_pos),(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])) / norm(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])^3
        else 
            return 0.0
        end
    end

    term_two = 0.0
    for particle in 1:path.n_particles
        centroid_pos = [sum(path.beads[bead,particle,dimension] for bead in 1:path.n_beads) for dimension in 1:path.n_dimensions]
        centroid_pos /= path.n_beads
        term_two += sum(get_term_two(particle, bead, other_bead, centroid_pos) for bead in 1:path.n_beads, other_bead in 1:path.n_beads)
    end

    return term_one - (t2_prefactor * term_two)
end

#multi threaded
function kinetic_energy(path::Path, potential::FrohlichPotential, estimator::Virial_EstimatorX)
    term_one = path.n_dimensions * path.n_particles / (2 * path.τ * path.n_beads)
    #term prefactor
    t2_prefactor = path.τ / (2 * path.n_beads)

    #beta conversion g_factor
    #β_conversion_factor = 47.81
    β_conversion_factor = 1

    β_reduced = path.τ * path.n_beads * β_conversion_factor


    function get_term_two(particle, bead, other_bead, centroid_pos)
        if bead != other_bead
            g_factor = 0.5 * potential.α * sqrt(potential.ħ/(2*path.m*potential.ω)) * cosh(β_reduced  * (abs(bead - other_bead)/path.n_beads - 0.5)) * csch(β_reduced * 0.5)
            return g_factor * dot((path.beads[bead,particle,:] - centroid_pos),(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])) / norm(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])^3
        else 
            return 0.0
        end
    end

    term_two = 0.0
    for particle in 1:path.n_particles
        centroid_pos = [ThreadsX.sum(path.beads[bead,particle,dimension] for bead in 1:path.n_beads) for dimension in 1:path.n_dimensions]
        centroid_pos /= path.n_beads
        term_two += ThreadsX.sum(get_term_two(particle, bead, other_bead, centroid_pos) for bead in 1:path.n_beads, other_bead in 1:path.n_beads)
    end

    return term_one - (t2_prefactor * term_two)
end


#Potential energy estimators -------------------------------------------


function potential_energy(path::Path, potential::OneBodyPotential, estimator::Estimator)
    return ThreadsX.sum(one_body_potential(potential, path, bead, particle)/path.n_beads for bead in 1:path.n_beads, particle in 1:path.n_particles)
end



function potential_energy(path::Path, potential::ConstantPotential, estimator::Thermodynamic_Estimator)
    return potential.V
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
    return kinetic_energy(path, potential, estimator) + potential_energy(path, potential, estimator)
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

# Position ------------------------------------------------------------------------


# Position
function Mean_Position(path::Path, potential::Potential, estimator::Estimator)
    positions = [mean(path.beads[:,particle,:]) for particle in path.n_particles]
    return positions
end

function Position(path::Path, potential::Potential, estimator::Estimator)
    return path.beads
end


