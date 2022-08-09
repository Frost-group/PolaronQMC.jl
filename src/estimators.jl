# estimators.jl
using Statistics

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
    term_one = path.n_dimensions * path.n_particles / (2*path.τ)
	prefactor = 1.0 / (4.0 * path.λ * path.τ^2)
	for bead in 1:path.n_beads, particle in 1:path.n_particles
		kinetic_energy += norm(path.beads[bead, particle, :] - path.beads[bead-1, particle, :])^2
	end
    return term_one - (kinetic_energy * prefactor / path.n_beads)
end

#=
function kinetic_energy(path::Path, estimator::Virial_Estimator)
    kinetic_energy = 0.0
    term_one = (path.n_dimensions * path.n_particles) / (2 * estimator.L * path.τ)
    prefactor = 1.0 / (4.0 * estimator.L * path.τ^2 * path.λ)
    for bead in 1:path.n_beads, particle in 1:path.n_particles
        link_term = dot(path.beads[bead+estimator.L, particle, :] - path.beads[bead, particle, :], path.beads[bead+1, particle, :] - path.beads[bead, particle, :])
        kinetic_energy += (term_one - prefactor * link_term)
    end
    return kinetic_energy / path.n_beads
end
=#

#=

function kinetic_energy(path::Path, estimator::Virial_Estimator)
    kinetic_energy = 0.0


    #term 1 in virial estimator for KE
    term_one = (path.n_dimensions * path.n_particles) / (2 * estimator.L * path.τ) 
    prefactor = 1.0 / (4.0 * estimator.L * path.τ^2 * path.λ)
    for bead in 1:path.n_beads, particle in 1:path.n_particles

        #term 2 in virial estimator for KE
            link_term = dot(path.beads[bead+estimator.L, particle, :] - path.beads[bead, particle, :], path.beads[bead+1, particle, :] - path.beads[bead, particle, :])
            term_two = prefactor*link_term
        kinetic_energy += term_two
    end
    return term_one - (kinetic_energy / path.n_beads)

end
=#

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
        Δ = zeros(path.n_dimensions)

        for j in -1*estimator.L+1:estimator.L-1
            Δ += (path.beads[bead, particle, :] - path.beads[bead+j, particle, :])
        end
        Δ ./= 2*estimator.L

        generalised_force = -1 / path.τ * path.m * potential.ω^2 * path.beads[bead, particle, :]
        potential_energy += -0.5*dot(generalised_force,Δ) + one_body_potential(potential, path, bead, particle)
    end
    return potential_energy / path.n_beads
end

#=
function potential_energy(path::Path, potential::FrohlichPotential, estimator::Virial_Estimator) #needs fixing
    prefactor_1 = potential.α * potential.ω^1.5 / (2*sqrt(2*path.m)*sinh(potential.β/2))
    βn = potential.β/path.n_beads
    potential_energy = 0.0
    
    for bead in 1:path.n_beads, particle in 1:path.n_particles
        Δ = zeros(path.n_dimensions)

        for j in -1*estimator.L+1:estimator.L-1
            Δ += path.beads[bead, particle, :] - path.beads[bead+j, particle, :]
        end
        Δ ./= 2*estimator.L


        #finding the derivative of the potential
        deriv_potential = zeros(path.n_dimensions) #derivative of potential wrt bead
        for other_bead in 1:path.n_beads
            if other_bead != bead
                prefactor = prefactor_1 * cosh(βn * abs(bead - other_bead) - potential.β/2)
                deriv_potential += prefactor .* (path.beads[bead,particle,:] - path.beads[other_bead,particle,:]) / norm(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])^3
            end
        end
        deriv_potential *= -1*βn


        #piecing all together
        generalised_force = -1 / path.τ * deriv_potential
        potential_energy += -0.5 * dot(generalised_force,Δ) + one_body_potential(potential, path, bead, particle)
    end
    return potential_energy / path.n_beads
end
=#



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



function Energy(path::Path, potential::Potential, estimator::Union{Thermodynamic_Estimator,Simple_Estimator})
    return kinetic_energy(path, estimator) + potential_energy(path, potential, estimator)
end


function Energy(path::Path, potential::FrohlichPotential, estimator::Virial_Estimator)
    term_one = path.n_dimensions * path.n_particles / (2 * path.τ * path.n_beads)
    #term prefactor
    t2_prefactor = path.τ / (2 * path.n_beads)

    #beta conversion g_factor
    β_conversion_factor = 47.81
    β_reduced = path.τ * path.n_beads * β_conversion_factor

    #needs to change to get working for multiple particles
    term_two = 0.0
    potential_energy = 0.0

    for particle in 1:path.n_particles
        centroid_pos = zeros(3)
        for bead in 1:path.n_beads
            centroid_pos += path.beads[bead,particle,:]
        end
        centroid_pos /= path.n_beads
    




        for bead in 1:path.n_beads
            potential_energy += one_body_potential(potential, path, bead, particle)

            for other_bead in 1:path.n_beads
                if other_bead != bead

                g_factor = 0.5 * potential.α * sqrt(potential.ħ/(2*path.m*potential.ω)) * cosh(β_reduced  * (abs(bead - other_bead)/path.n_beads - 0.5)) * csch(β_reduced * 0.5)
                term_two += g_factor * dot((path.beads[bead,particle,:] - centroid_pos),(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])) / norm(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])^3
                end
            end

            
        end
    end

    #can extend to sum up energy of all particles 

    return term_one - (t2_prefactor * term_two) + (potential_energy / path.n_beads)
    
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
    return path.beads[:,1,1]
end


