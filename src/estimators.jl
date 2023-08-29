# estimators.jl
using Statistics
using ThreadsX
using SpecialFunctions



#-------------Kinetic energy estimators-----------------------
function kineticEnergy(path::Path, potential::Potential, estimator::Union{SimpleEstimator, SimpleVirialEstimator}) #thermal dynamic estimator from ceperly paper
	kinetic_energy = 0.0

    # Add-up energy for each particle and each bead
	for bead in 1:path.n_beads, particle in 1:path.n_particles
		kinetic_energy += 0.5 * path.m[particle] * norm(path.beads[bead, particle] - path.beads[mod1(bead-1, path.n_beads), particle])^2 / path.τ^2
	end

	return kinetic_energy / path.n_beads
end

function kineticEnergy(path::Path, potential::Potential, estimator::ThermodynamicEstimator, store_diff::Vector{Float64}) #thermal dynamic estimator from ceperly paper
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


function kineticEnergy(path::Path, potential::Union{HarmonicPotential, HarmonicInteractionPotential}, estimator::VirialEstimator, store_diff::Vector{Float64})
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


function kineticEnergy(path::Path, potential::FrohlichPotential, estimator::VirialEstimator, store_diff::Vector{Float64})
    """
    Virial estimator from Ceperly's book Interactiong Electrons
    For Frohlich Oscillator

    Relate the average kinetic energy with potential energy -> Hence used centroid position for EACH PARTICLE
    """
    n_beads = path.n_beads;
    β = path.τ * n_beads;

    term_one = (path.n_dimensions * path.n_particles) / (2 * path.τ * n_beads) # same as harmonic
    final_term_two = 0.0;
    #Defining constants

    function getTermTwo(particle, bead, other_bead, centroid_pos, ω)
        if bead != other_bead
            g_factor = cosh(potential.ħ * ω * β * (abs(bead-other_bead)/n_beads - 0.5))
            for i in 1:path.n_dimensions
                store_diff[i] = path.beads[bead,particle,i] - path.beads[other_bead,particle,i]
            end
            #return g_factor * dot((path.beads[bead,particle,:] - centroid_pos),(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])) / norm(path.beads[bead,particle,:] - path.beads[other_bead,particle,:])^3
            return g_factor * dot((@views path.beads[bead,particle,:] - centroid_pos), store_diff) / norm(store_diff)^3
        else 
            return 0.0
        end
    end

    for ω in potential.ω
        term_two = 0.0
        #t2_prefactor = path.τ / (2 * (path.n_beads-1)) * 0.5 * potential.α * (potential.ħ * potential.ω)^(3/2) * sqrt(1/2/path.m) * csch(potential.ħ * potential.ω * β / 2)
        t2_prefactor = path.τ/(2 * (n_beads)) * 0.5 * potential.α * (potential.ħ * ω)^(3/2) * sqrt(1/2/path.m) * csch(potential.ħ * ω * β / 2)
    
        # F = -dV/dr ∝ r/(r(τ)-r(τ'))^3 (Classical force formula from the pseudo-potential)

        for particle in 1:path.n_particles
            centroid_pos = [sum(path.beads[bead,particle,dimension] for bead in 1:n_beads) for dimension in 1:path.n_dimensions] / n_beads
            #centroid_pos /= path.n_beads
            term_two += sum(getTermTwo(particle, bead, other_bead, centroid_pos, ω) for bead in 1:n_beads, other_bead in 1:n_beads)
        end

        final_term_two += term_two * t2_prefactor;
    end

    return term_one + final_term_two
    #return term_one + (t2_prefactor * term_two) # -1 (from eqn) * -1 (frm dV/dr) * -1 (force formula) * (-1) from potential [updated]
    #return (t2_prefactor * term_two) # -1 (from eqn) * -1 (frm dV/dr) * -1 (force formula) * (-1) from potential [updated]
end

function kineticEnergy(path::DiscretePath, potential::HolsteinPotential, estimator::VirialEstimator)
    kinetic_energy = 0.0
    for particle in 1:path.n_particles
        for i in 1:path.n_beads
            for k in 1:path.n_dimensions
                kinetic_energy += (besseli(path.beads[i, particle, k] - path.beads[mod1(i+1, path.n_beads), particle, k] - 1, 2 * path.τ * potential.J) + besseli(path.beads[i, particle, k] - 
                                path.beads[mod1(i+1, path.n_beads), particle, k] + 1, 2 * path.τ * potential.J))/(besseli(path.beads[i, particle, k] - path.beads[mod1(i+1, path.n_beads), particle, k], 2 * path.τ * potential.J))
        
            end
        end
        
        #=
        kinetic_energy += (besseli(path.beads[i] - path.beads[mod1(i+1, path.n_beads)] - 1, 2 * path.τ * potential.J) + besseli(path.beads[i] - 
                            path.beads[mod1(i+1, path.n_beads)] + 1, 2 * path.τ * potential.J))/(besseli(path.beads[i] - path.beads[mod1(i+1, path.n_beads)], 2 * path.τ * potential.J))
    
        =#
    end
    return -1 * kinetic_energy/path.n_beads
end


#-------------Potential energy estimators---------------------
#=
function potentialEnergy(path::Path, potential::FrohlichPotential, estimator::Estimator, store_diff::Vector{Float64}, prop_Matrix::Array{Float64})
    potential_energy = 0.0
    n_beads, l_ω = path.n_beads, length(potential.ω);
    β, α = path.τ * n_beads, potential.α;
    #term_factor = -0.5 * α * (ħω)^(3/2) * sqrt(1/2/path.m) * csch(ħω * β / 2);

    # Calculates the double integral component with individual modes

    for particle in 1:path.n_particles
        for bead in 1:n_beads
            for other_bead in 1:n_beads
                if other_bead != bead # Discounting self-contributions
                    for i in 1:path.n_dimensions
                        store_diff[i] = path.beads[mod1(bead, n_beads), particle, i] - path.beads[mod1(other_bead, n_beads), particle, i]
                    end

                    for i in 1:l_ω
                        #=
                        potential_energy += -0.5 * α * (potential.ħ * potential.ω[i])^(3/2) * sqrt(1/2/path.m) /norm(store_diff) *
                        (path.τ^2 * sinh(potential.ħ * β * potential.ω[i] * (abs(other_bead-bead)/path.n_beads-0.5)) * csch(potential.ħ * potential.ω[i] * β / 2) * potential.ħ * potential.ω[i] * path.n_beads * (abs(other_bead-bead)/path.n_beads-0.5)
                        - path.τ^2 * cosh(potential.ħ * potential.ω[i] * β / 2) /sinh(potential.ħ * potential.ω[i] * β / 2)^2 * prop_Matrix[max(bead, other_bead), min(bead, other_bead), i] * 0.5 * potential.ħ * potential.ω[i] * path.n_beads
                        + 2 * path.τ * prop_Matrix[max(bead, other_bead), min(bead, other_bead), i] * csch(potential.ħ * potential.ω[i] * β / 2))
                        =#
                        potential_energy += -0.5 * α * (potential.ħ * potential.ω[i])^(3/2) * sqrt(1/2/path.m) /norm(store_diff) * 2 * path.τ * prop_Matrix[max(bead, other_bead), min(bead, other_bead), i] * csch(potential.ħ * potential.ω[i] * β / 2)

                    end
                    
                end
            end
        end
    end

    return potential_energy/path.n_beads

end
=#

function potentialEnergy(path::Path, potential::OneBodyPotential, estimator::Estimator, store_diff::Vector{Float64}, prop_Matrix::Array{Float64})
    return sum(oneBodyPotential(potential, path, bead, particle, store_diff, prop_Matrix)/path.n_beads for bead in 1:path.n_beads, particle in 1:path.n_particles)
end

function potentialEnergy(path::DiscretePath, potential::HolsteinPotential, estimator::Estimator, DF_l::Array{Float64})
    potential_energy = 0.0
    for particle in 1:path.n_particles
        for i in 1:path.n_beads
            for j in 1:path.n_beads
                # Delta function to check whether the 2 beads are on the same sites
                if δ(@view(path.beads[i, particle, :]), @view(path.beads[j, particle, :]))
                    potential_energy += DF_l[abs(i-j)+1]
                end

            end
        end
    end

    
    return -potential_energy/path.n_beads
    #return sum(oneBodyPotential(potential, path, bead, particle, store_diff, prop_Matrix)/path.n_beads for bead in 1:path.n_beads, particle in 1:path.n_particles)
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


function energy(path::Path, potential::Potential, estimator::Estimator, store_diff::Vector{Float64}, prop_Matrix::Array{Float64})
    """
    return kineticEnergy(path, potential, estimator) + potentialEnergy(path, potential, estimator)
    """
    #Printing out energy for quick inspection -> Sometimes the values diverges to -Inf and can stop at early stage
    KE::Float64 = kineticEnergy(path, potential, estimator, store_diff)
    PE::Float64 = potentialEnergy(path, potential, estimator, store_diff, prop_Matrix)
    println("KE:", trunc(KE, digits=2), " ", "PE:", trunc(PE, digits=2))
    if typeof(potential) == FrohlichPotential
        return KE + PE - 1.5 * 1/(exp(path.τ * path.n_beads) - 1)
    else
        return KE + PE
    end
end

function energy(path::DiscretePath, potential::Potential, estimator::Estimator, F_l::Array{Float64})
    """
    return kineticEnergy(path, potential, estimator) + potentialEnergy(path, potential, estimator)
    """
    #Printing out energy for quick inspection -> Sometimes the values diverges to -Inf and can stop at early stage
    KE::Float64 = kineticEnergy(path, potential, estimator)
    PE::Float64 = potentialEnergy(path, potential, estimator, F_l)
    #println("KE:", trunc(KE, digits=2), " ", "PE:", trunc(PE, digits=2))
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


