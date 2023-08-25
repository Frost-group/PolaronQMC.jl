# actions.jl
# Ways of calculating action
using SpecialFunctions

# usually the total action is not used because we calcalate changes instead of the full action. The formula is put here for completeness
function totalAction(path::Path, bead_one::Int, bead_two::Int, particle::Int, potential::OneBodyPotential, regime::Regime)
    return kinecticAction(path, bead_one, bead_two, particle, regime) + potential_action(path, bead_one, particle, potential, regime)
end

"""
    Calculate Kinetic Actions; same across all regime, just slightly different formalism

    KineticAction(path::Path, bead_one::Int, bead_two::Int, particle::Int)

    Kinetic action for a link between two beads `bead_one` and `bead_two` for a particle indexed by an integer `particle`. Usually, these are neighbouring beads.

    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead_one::Int`: first bead in the link.
    - `bead_two::Int`: second bead in the link.
    - `particle::Int`: select a specific particle indexed by this integer.

    See also [`Path`](@ref). 
"""

function kineticAction(path::Path, bead_one::Int, bead_two::Int, particle::Int, regime::Regime, potential::Potential, store_diff::Vector{Float64})
    """
    Allow multiple particles with different mass
    """
    kinetic_action::Float64 = 0.5 * path.m[particle] * norm(path.beads[mod1(bead_two, path.n_beads), particle, :] - path.beads[mod1(bead_one, path.n_beads), particle, :])^2 / path.τ
    return kinetic_action
end

#Primitive method (based off Ceperly paper)
function kineticAction(path::Path, bead_one::Int64, bead_two::Int64, particle::Int64, regime::PrimitiveRegime, potential::Potential, store_diff::Vector{Float64})
    """
    Allow multiple particles with different mass/λ
    """
    #kinetic_action::Float64 = norm(path.beads[mod1(bead_two, path.n_beads), particle, :] - path.beads[mod1(bead_one, path.n_beads), particle, :])^2 / (4 * path.λ[particle] * path.τ)
    for i in 1:path.n_dimensions
        store_diff[i] = path.beads[mod1(bead_two, path.n_beads), particle, i] - path.beads[mod1(bead_one, path.n_beads), particle, i]
    end
    kinetic_action = norm(store_diff)^2 / path.K_factor
    #@timeit tmr "Kinetic" kinetic_action = norm(path.beads[mod1(bead_two, path.n_beads), particle, :] - path.beads[mod1(bead_one, path.n_beads), particle, :])^2 / path.K_factor
    return kinetic_action
end

function kineticAction(path::DiscretePath, bead_one::Int64, bead_two::Int64, particle::Int64, potential::HolsteinPotential)
    # Thermodynamic limit N -> ∞
    kinetic_action = 1.0;
    for dimension in 1:path.n_dimensions
        kinetic_action *= besseli(path.beads[bead_one, particle, dimension] - path.beads[bead_two, particle, dimension], 2 * path.τ * potential.J)
    end
    return kinetic_action
end




"""
Calculate Potential Action
"""

function potentialAction(path::Path, bead::Int, particle::Int, potential::OneBodyPotential, regime::SimpleRegime)
    """
    Calculate potential action {For One-Body Potential, Simple Regime}

    potentialAction(path::Path, bead_one::Int, bead_two::Int, particle::Int)

    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead::Int`: bead of interest.
    - `particle::Int`: select a specific particle indexed by this integer.
    - `potentials::OneBodyPotential`: any one-body potential.

    See also [`Path`](@ref). 
    """
    return oneBodyPotential(potential, path, bead, particle) * path.τ
end

function potentialAction(path::Path, bead::Int, particle::Int, potential::OneBodyPotential, regime::PrimitiveRegime, store_diff::Vector{Float64}, prop_Matrix::Array{Float64})

    """
    Calculate potential action {For One-body Potential, Primitive Regime}

    potential_action(path::Path, bead::Int, particle::Int, potential::OneBodyPotential)

    Potential action for a bead indexed by an integer `bead` and a particle indexed by an integer `particle`.
    This is the potential action for an external one-body potential.

    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead::Int`: select a specific bead indexed by an integer.
    - `particle::Int`: select a specific particle indexed by this integer.
    - `potentials::OneBodyPotential`: a one-body potential type. Acts like an external potential on the system.

    See also [`Path`](@ref), [`OneBodyPotential`](@ref). 
    """
    
    #@timeit tmr "PotentialF" 
    if typeof(potential) == FrohlichPotential
        # Extra factor of 2
        return 2 * path.τ * oneBodyPotential(potential, path, bead, particle, store_diff, prop_Matrix)
    end

    #@timeit tmr "Potential"
    return path.τ * oneBodyPotential(potential, path, bead, particle, store_diff, prop_Matrix)
end

function potentialAction(path::DiscretePath, bead::Int, particle::Int, potential::HolsteinPotential, F_l::Array{Float64})
    return oneBodyPotential(potential, path, bead, particle, F_l)
end

function potentialAction(path::Path, bead::Int, particle::Int, potential::ConstantPotential, regime::Regime)

    """
    Calculate potential action {For Constant Potential, Any Regime}

    potential_action(path::Path, bead::Int, particle::Int, potential::ConstantPotential)

    Potential action for a bead indexed by an integer `bead` and a particle indexed by an integer `particle`.
    This is the potential action for an external constant potential.

    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead::Int`: select a specific bead indexed by an integer.
    - `particle::Int`: select a specific particle indexed by this integer.
    - `potentials::ConstantPotential`: a constant potential type. Typically just zero.

    See also [`Path`](@ref), [`ConstantPotential`](@ref). 
    """
    # Potential energy is V, but action is τV 
    return path.τ * potential.V
end

function potentialAction(path::Path, bead::Int, particle::Int, potential::TwoBodyPotential, regime::Regime)

    """
    Calculate potential action {For Two-body Potential, Any Regime}

    potentialAction(path::Path, bead::Int, particle::Int, potential::TwoBodyPotential)

    Potential action for a bead indexed by an integer `bead` and a particle indexed by an integer `particle`.
    This is the potential action for an internal two-body potential between the specified bead on the specified particle and the same bead on all other particles. 

    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead::Int`: select a specific bead indexed by an integer.
    - `particle::Int`: select a specific particle indexed by this integer.
    - `potentials::TwoBodyPotential`: a two-body potential type. Acts like interactions between pairs of particles in the system.

    See also [`Path`](@ref), [`TwoBodyPotential`](@ref). 
    """

    potential_action = sum(twoBodyPotential(potential, path, bead, particle, other_particle) for other_particle in 1:path.n_particles if particle != other_particle)
    return path.τ * potential_action 
end

function potentialAction(path::Path, bead::Int, particle::Int, potential::FrohlichPotential, regime::LBRegime)
    """
    Soecifically written Li-Broughton regime for more accurate potential calculation

    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead::Int`: select a specific bead indexed by an integer.
    - `particle::Int`: select a specific particle indexed by this integer.
    - `potentials::OneBodyPotential`: a one-body potential type. Acts like an external potential on the system.

    Ref TESI paper
    """
    # Define constants, avoid repeated calling
    β = path.τ * path.n_beads
    t2_prefactor = path.τ^3 / (24 * path.m)
    V2_prefactor = (0.5 * potential.α * (potential.ħ * potential.ω)^(3/2) * sqrt(1/2/path.m) * csch(potential.ħ * potential.ω * β / 2))^2
    ħω = potential.ω * potential.ħ
    # F = -dV/dr ∝ r/(r(τ)-r(τ'))^3
    
    inner_integral = 0.0
    for other_bead in 1:path.n_beads+1
        if mod1(bead, path.n_beads) != mod1(other_bead, path.n_beads)
            g_factor = (cosh(ħω * β * (abs(bead-other_bead)/(path.n_beads) - 0.5)))^2
            inner_integral += g_factor / norm(path.beads[mod1(bead, path.n_beads),particle,:] - path.beads[mod1(other_bead, path.n_beads),particle,:])^4
        end
    end
    return 2 * path.τ * oneBodyPotential(potential, path, bead, particle) + 2 * inner_integral * path.τ^2 * t2_prefactor * V2_prefactor
end

function potentialAction(path::Path, bead::Int, beadrange::Union{Vector{Int64}, StepRange{Int64, Int64}, UnitRange{Int64}}, particle::Int, potential::FrohlichPotential, regime::PrimitiveRegime)
    """
    Specifically written for Frohlich in bisectMover due to the complexity induced by double time integrals

    bisectPotentialAction(path::Path, bead::Int, beadrange::Union{Vector{Int64}, StepRange{Int64, Int64}, UnitRange{Int64}}, particle::Int, potential::FrohlichPotential, regime::PrimitiveRegime)
    
    Arguments
    - `path::Path`: collection of all particle imaginary-time paths.
    - `bead::Int`: bead of interest.
    - `beadrange::Union{Vector{Int64}, StepRange{Int64, Int64}, UnitRange{Int64}}`: bead within the same segment to avoid double counting contributions
    - `particle::Int`: select a specific particle indexed by this integer.
    - `potentials::OneBodyPotential`: any one-body potential.

    """

    # Defining constants (DON'T have to call path.variable iteratively)
    β = path.τ * path.n_beads; m = path.m; ħω = potential.ħ * potential.ω; α = potential.α;
    n_beads = path.n_beads
    term_factor = -0.5 * α * (ħω)^(3/2) * sqrt(1/2/m) * csch(ħω * β / 2)
    bead = mod1(bead, n_beads)

    # Integrals contribution summation
    inner_integral = 0.0
    for other_bead in 1:n_beads
        if other_bead != bead
            if other_bead in beadrange
                inner_integral += 0.5 * cosh(ħω * β * (abs(bead-other_bead)/n_beads - 0.5)) / norm(path.beads[mod1(bead, n_beads), particle, :] - path.beads[mod1(other_bead, n_beads), particle, :])
            else
                inner_integral += cosh(ħω * β * (abs(bead-other_bead)/n_beads - 0.5)) / norm(path.beads[mod1(bead, n_beads), particle, :] - path.beads[mod1(other_bead, n_beads), particle, :])
            end
        end
    end
    return 2 * (path.τ)^2 * inner_integral * term_factor # Note that this path.τ multiplication refer to dτ'
end
