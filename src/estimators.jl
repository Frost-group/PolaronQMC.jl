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
     #term_one = (path.n_dimensions * path.n_particles) / (path.τ * path.n_beads) 
 
     #term 2 prefactor
     t2_prefactor = (path.m * potential.ω^2) / (2 * path.n_beads)
     term_two = 0.0
     for particle in 1:path.n_particles
         centroid_pos = [ThreadsX.sum(path.beads[bead,particle,dimension] for bead in 1:path.n_beads) for dimension in 1:path.n_dimensions]
         centroid_pos /= path.n_beads
 

        term_two = ThreadsX.sum(dot(path.beads[bead, particle, :] - centroid_pos, path.beads[bead, particle, :]) for bead in 1:path.n_beads)

     end
 
     return term_one + t2_prefactor * term_two
 
 end


function kinetic_energy(path::Path, potential::FrohlichPotential, estimator::Virial_Estimator)
    β = path.τ * path.n_beads

    #term prefactor
    #term_one = 0.0 #Not sure why is 0
    term_one = (path.n_dimensions * path.n_particles) / (2 * path.τ * path.n_beads) # same as harmonic
    t2_prefactor = path.τ / (2 * path.n_beads) * 0.5 * potential.α * (potential.ħ * potential.ω)^(3/2) * sqrt(1/2/path.m) * csch(potential.ħ * potential.ω * β / 2)
    ħω = potential.ω * potential.ħ
    # F = -dV/dr ∝ r/(r(τ)-r(τ'))^3
    
    function get_term_two(particle, bead, other_bead, centroid_pos)
        if bead != other_bead
            #g_factor = -0.5 * potential.α * (potential.ħ * potential.ω)^3/2 * sqrt(2*path.m) * cosh(potential.ω*β * (abs(bead-other_bead)/path.n_beads - 0.5 * potential.ħ)) * csch(potential.ħ * potential.ω * β / 2)
            #g_factor = 0.5 * potential.α * (potential.ħ * potential.ω)^3/2 / pi * sqrt(1/2/path.m) * cosh(potential.ω * β * potential.ħ * (abs(bead-other_bead)/path.n_beads - 0.5))* csch(potential.ħ * potential.ω * β / 2)
            #g_factor = 0.5 * potential.α * (potential.ħ * potential.ω)^(3/2) * sqrt(1/2/path.m) * cosh(potential.ω * β * potential.ħ * (abs(bead-other_bead)/path.n_beads - 0.5)) * csch(potential.ħ * potential.ω * β / 2)
            g_factor = cosh(ħω * β * (abs(bead-other_bead)/path.n_beads - 0.5))
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

    #return term_one - (t2_prefactor * term_two) #original one from George
    return term_one + (t2_prefactor * term_two) # -1 (from eqn) * -1 (frm dV/dr) * -1 (force formula) * (-1) from potential [updated]
    #return term_one - (t2_prefactor * term_two) + (1.5 * potential.ħ * potential.ω * coth(potential.ω * potential.ħ * β)) #Phonon Kinetic Energy
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
    KE = kinetic_energy(path, potential, estimator)
    PE = potential_energy(path, potential, estimator)
    println("kinetic is:", KE)
    println("Potential is:", PE)
    return KE + PE
    #return kinetic_energy(path, potential, estimator) + potential_energy(path, potential, estimator)
end


# Correlation ---------------------------------------------------------------------

function Correlation(path::Path, potential::Potential, estimator::Estimator)
    correlation = Vector{Float64}(undef, path.n_beads-1)
    #=
    for Δτ in 1:path.n_beads, bead_one in 1:path.n_beads
        bead_two = mod1(bead_one + Δτ, path.n_beads)
        correlation[Δτ] += dot(path.beads[bead_one, :, :], path.beads[bead_two, :, :])
    end
    =#
    for Δτ in 1:(path.n_beads-1)
        
        for bead_one in 1:path.n_beads
            bead_two = bead_one + Δτ
            if bead_two <= path.n_beads
                correlation[Δτ] += dot(path.beads[bead_one, :, :], path.beads[bead_two, :, :])
            end
        end
        correlation[Δτ] /= (path.n_beads - Δτ)
    end				
    return [correlation]
end

function autoCorrelation(path::Path, potential::Potential, estimator::Estimator)

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


