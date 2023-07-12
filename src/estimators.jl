# estimators.jl
using Statistics
using ThreadsX


#Kinectic energy estimators ---------------------------------

function kineticEnergy(path::Path, potential::Potential, estimator::Union{SimpleEstimator, SimpleVirialEstimator}) #thermal dynamic estimator from ceperly paper
	kinetic_energy = 0.0

    # Add-up energy for each particle and each bead
	for bead in 1:path.n_beads, particle in 1:path.n_particles
		kinetic_energy += 0.5 * path.m[particle] * norm(path.beads[bead, particle] - path.beads[mod1(bead-1, path.n_beads), particle])^2 / path.τ^2
	end

	return kinetic_energy / path.n_beads
end



function kineticEnergy(path::Path, potential::Potential, estimator::ThermodynamicEstimator) #thermal dynamic estimator from ceperly paper
    """
    Thermodynamic estimator from Ceperly's book Interactiong Electrons
    """
    kinetic_energy = 0.0 # For the bead kinetic energy
    term_one = path.n_dimensions * path.n_particles / (2*path.τ) # Classical energy term, N/2*(kT) (constant)
	
    prefactor = 1.0 / (4.0 * path.λ * path.τ^2)
	for bead in 1:path.n_beads, particle in 1:path.n_particles
		kinetic_energy += norm(path.beads[bead, particle, :] - path.beads[mod1(bead-1, path.n_beads), particle, :])^2
	end
    return term_one - (kinetic_energy * prefactor / path.n_beads)
end


function kineticEnergy(path::Path, potential::Union{HarmonicPotential, HarmonicInteractionPotential}, estimator::VirialEstimator)
    """
    Virial estimator from Ceperly's book Interactiong Electrons
    For Harmonic Oscillator

    Relate the average kinetic energy with potential energy -> Hence used centroid position for EACH PARTICLE
    """
    term_one = (path.n_dimensions * path.n_particles) / (2 * path.τ * path.n_beads) # Classical energy term (constant)
 
    #term_2 prefactor
    t2_prefactor = (path.m * potential.ω^2) / (2 * path.n_beads)
    term_two = 0.0
    for particle in 1:path.n_particles
        centroid_pos = [ThreadsX.sum(path.beads[mod1(bead, path.n_beads),particle,dimension] for bead in 1:path.n_beads) for dimension in 1:path.n_dimensions]
        centroid_pos /= path.n_beads

        term_two += ThreadsX.sum(dot(path.beads[mod1(bead, path.n_beads), particle, :] - centroid_pos, path.beads[mod1(bead, path.n_beads), particle, :]) for bead in 1:path.n_beads)
    end
 
    return term_one + t2_prefactor * term_two
end


function kineticEnergy(path::Path, potential::FrohlichPotential, estimator::VirialEstimator)
    """
    Virial estimator from Ceperly's book Interactiong Electrons
    For Frohlich Oscillator

    Relate the average kinetic energy with potential energy -> Hence used centroid position for EACH PARTICLE
    """

    #Defining constants
    β = path.τ * path.n_beads
    term_one = (path.n_dimensions * path.n_particles) / (2 * path.τ * path.n_beads) # same as harmonic
    t2_prefactor = path.τ / (2 * path.n_beads) * 0.5 * potential.α * (potential.ħ * potential.ω)^(3/2) * sqrt(1/2/path.m) * csch(potential.ħ * potential.ω * β / 2)
    ħω = potential.ω * potential.ħ
    # F = -dV/dr ∝ r/(r(τ)-r(τ'))^3 (Classical force formula from the pseudo-potential)
    
    function getTermTwo(particle, bead, other_bead, centroid_pos)
        if bead != other_bead
            g_factor = cosh(ħω * β * (abs(bead-other_bead)/path.n_beads - 0.5))
            #A = 1/(norm(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])^3)
            #if A != Inf
            return g_factor * dot((path.beads[bead,particle,:] - centroid_pos),(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])) / norm(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])^3
                #return g_factor * dot((path.beads[bead,particle,:] - centroid_pos),(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])) * A
            #end
        else 
            return 0.0
        end
    end

    term_two = 0.0
    for particle in 1:path.n_particles
        centroid_pos = [sum(path.beads[bead,particle,dimension] for bead in 1:path.n_beads) for dimension in 1:path.n_dimensions]
        centroid_pos /= path.n_beads
        term_two += sum(getTermTwo(particle, bead, other_bead, centroid_pos) for bead in 1:path.n_beads, other_bead in 1:path.n_beads)
    end

    return term_one + (t2_prefactor * term_two) # -1 (from eqn) * -1 (frm dV/dr) * -1 (force formula) * (-1) from potential [updated]
end

#Potential energy estimators -------------------------------------------

function potentialEnergy(path::Path, potential::OneBodyPotential, estimator::Estimator)
    return sum(oneBodyPotential(potential, path, bead, particle)/(path.n_beads) for bead in 1:path.n_beads, particle in 1:path.n_particles)
end


function potentialEnergy(path::Path, potential::ConstantPotential, estimator::ThermodynamicEstimator)
    return potential.V
end


function potentialEnergy(path::Path, potential::TwoBodyPotential, estimator::ThermodynamicEstimator)
    potential_energy = 0.0
    for bead in 1:path.n_beads, particle_one in 1:path.n_particles, particle_two in 1:path.n_particles
        if particle_one != particle_two
            potential_energy += twoBodyPotential(potential, path, bead, particle_one, particle_two)
        end
    end
    return potential_energy / path.n_beads
end


# Energy ---------------------------------------------------------


function energy(path::Path, potential::Potential, estimator::Estimator)
    #return kineticEnergy(path, potential, estimator) + potentialEnergy(path, potential, estimator)
    #Printing out energy for quick inspection -> Sometimes the values diverges to -Inf and can stop at early stage
    KE = kineticEnergy(path, potential, estimator)
    PE = potentialEnergy(path, potential, estimator)
    println("kinetic is:", KE, " ", "Potential is:", PE)
    #println("Potential is:", PE)
    return KE + PE
end


# Correlation ---------------------------------------------------------------------

function correlation(path::Path, potential::HarmonicPotential, estimator::Estimator)
    """
    Position Correlation function from Vvedensky's paper (user's guide for Path-integral Monte Carlo)
    G(Δτ) = <x(τ)x(τ+Δτ)> - <x(τ)><x(τ+Δτ)>
    It is used to determine whether the τ is too small or too large
    
    For Harmonic oscillator it is in special form as <x(τ)> = 0
        
    output
        correlation: has array of size n_beads-1 since the maximum Δτ difference is given by (n-1)*Δτ
    """
    
    correlation = zeros(path.n_beads-1)

    for Δτ in 1:(path.n_beads-1)
        for bead_one in 1:path.n_beads
            bead_two = bead_one + Δτ
            if bead_two > path.n_beads
                break
            end

            correlation[Δτ] += dot(path.beads[bead_one, :, :], path.beads[bead_two, :, :])
        end
        
        correlation[Δτ] /= (path.n_beads - Δτ)
    end				
    return correlation
end


function correlation(path::Path, potential::Potential, estimator::Estimator)
    """
    Position Correlation function from Vvedensky's paper (user's guide for Path-integral Monte Carlo)
    G(Δτ) = <x(τ)x(τ+Δτ)> - <x(τ)><x(τ+Δτ)>
    It is used to determine whether the τ is too small or too large
        
    output
        correlation: has array of size n_beads-1 since the maximum Δτ difference is given by (n-1)*Δτ
    """
    #Initialise correlation array
    correlation = zeros(path.n_beads-1)
    for Δτ in 1:(path.n_beads-1)
        avg_i = 0.0
        avg_iτ = 0.0
        for bead_one in 1:path.n_beads
            bead_two = bead_one + Δτ
            if bead_two > path.n_beads
                break
            end

            if avg_i == 0.0
                avg_i = path.beads[bead_one, :, :]
                avg_iτ = path.beads[bead_two, :, :]
            else
                avg_i += path.beads[bead_one, :, :]
                avg_iτ += path.beads[bead_two, :, :]
            end
            correlation[Δτ] += dot(path.beads[bead_one, :, :], path.beads[bead_two, :, :])
        end
        correlation[Δτ] /= (path.n_beads - Δτ)
        correlation[Δτ] -= dot(avg_i, avg_iτ) / ((path.n_beads - Δτ)^2)
    end				
    return correlation
end


function position(path::Path, potential::Potential, estimator::Estimator)
    return path.beads
end


