# types.jl

using StaticArrays


# Regime for the simuation to run in
abstract type Regime end


# Simple form of calculating action
struct SimpleRegime <: Regime  
    function SimpleRegime()
        new()
    end
end


# Calculating using the primitive approximation as per Ceperly paper
struct PrimitiveRegime <: Regime 
    function PrimitiveRegime()
        new()
    end
end

struct LBRegime <: Regime
    function LBRegime()
        new()
    end
end


mutable struct Path

    """
    Generic path mutable type 

    Attributes
        n_beads::Int: Number of beads from T and τ
        n_particles::Int: Number of particles
        n_dimensions::Int: Number of dimensions (1D, 2D, 3D)

        beads: Sized array that are used to store all beads positions 

        τ::Float64: size of each imaginary time slice
        m::Union{Float64, Vector{Float64}}: mass of the particles, can be a vector if we have multiparticles with different mass
        λ::Float64: simplify factor of ħ^2/2m

    """

	n_beads :: Int64
	n_particles :: Int64
    n_dimensions :: Int64

	beads :: SizedArray

	τ :: Float64
    m :: Union{Float64, Vector{Float64}} 
	λ :: Float64

	function Path(n_beads::Int64, n_particles::Int64, n_dimensions::Int64, τ::Float64; m = 1.0, λ = 0.5, start_range = 1.0)
        
        # Randomised initial bead position around 0. Use static array since we won't be modifying the number of beads in the simulation
        beads = @SArray randn(n_beads, n_particles, n_dimensions)
        beads *= start_range # If we want the beads to be more widespread

		new(n_beads, n_particles, n_dimensions, beads, τ, m, λ)
	end
end

# abstract type for how to move the beads in the simulation
abstract type Mover end


mutable struct SingleMover <: Mover
    """
    Displace one bead at a time
    """
    adjusters :: SizedArray

    function SingleMover(path::Path)
        particles = [SingleAdjuster(path.λ, path.τ) for i in 1:path.n_particles]
        adjusters = SVector(Tuple(particles))
        new(adjusters)
    end
end


mutable struct DisplaceMover <: Mover
    """
    Displace the whole chain at a time
    """
    adjusters :: SizedArray

    function DisplaceMover(path::Path)
        particles = [DisplaceAdjuster(path.λ, path.τ) for i in 1:path.n_particles]
        adjusters = SVector(Tuple(particles))
        new(adjusters)
    end
end


mutable struct BisectMover <: Mover
    """
    Displace according to mid-point strategies, layer by layer until whole segment moved
    """

    adjusters :: SizedArray

    function BisectMover(path::Path)
        particles = [BisectAdjuster(path.λ, path.τ) for i in 1:path.n_particles]
        adjusters = SVector(Tuple(particles))
        new(adjusters)
    end
end


# Adjuster type to container information about shift width and allow for its auto adjustment
abstract type Adjuster end


# Adjuster for the Single! move algorithm
mutable struct SingleAdjuster <: Adjuster

    """
    Adjuster for the Single! move algorithm
    """

    attempt_counter :: Int
    success_counter :: Int
    value :: Float64
    acceptance_rate :: Float64
    function SingleAdjuster(λ::Float64, τ::Float64)
        value = sqrt(4 * λ * τ) * 0.5
        new(0, 0, value, 0)
    end
end


#Adjuster for the Displace! move algorithm
mutable struct DisplaceAdjuster <: Adjuster

    """
    Adjuster for the Displace! move algorithm
    """

    attempt_counter :: Int
    success_counter :: Int
    value :: Float64
    acceptance_rate :: Float64
    function DisplaceAdjuster(λ::Float64, τ::Float64)
        value = 1
        new(0, 0, value, 0)
    end
end


#Adjuster for the Bisect! move algorithm
mutable struct BisectAdjuster <: Adjuster
    
    """
    Adjuster for the Bisect! move algorithm
    """

    attempt_counter :: Int
    success_counter :: Int
    value :: Float64
    acceptance_rate :: Float64

    function BisectAdjuster(λ::Float64, τ::Float64)
        value = sqrt(τ * λ)
        new(0, 0, value, 0)
    end
end


# Most abstract Potential type.
abstract type Potential end

# A potential that is independent on bodies.
abstract type NoBodyPotential <: Potential end

# For potentials that depend on one body.
abstract type OneBodyPotential <: Potential end

# For potentials that depend on two bodies.
abstract type TwoBodyPotential <: Potential end


struct ConstantPotential <: NoBodyPotential

    """
    A constant potential
    """

    V :: Float64
    function ConstantPotential(V::Float64)
        new(V)
    end
end


struct HarmonicPotential <: OneBodyPotential

    """
    A Harmonic potential for a single body.
    """

    ω :: Float64
    function HarmonicPotential(ω::Float64)
        new(ω)
    end
end


struct FrohlichPotential <: OneBodyPotential

    """
    Potential for Frohlich Polaron
    """

    α :: Float64
    ω :: Float64
    ħ :: Float64
    function FrohlichPotential(α::Float64, ω::Float64, ħ::Float64)
        new(α, ω, ħ)
    end
end


struct HarmonicInteractionPotential <: OneBodyPotential

    """
    A Harmonic potential for a single body.
    """

    ω :: Float64
    κ :: Float64

    function HarmonicInteractionPotential(ω::Float64, κ :: Float64)
        new(ω, κ)
    end
end


struct FrohlichInteractionPotential <: OneBodyPotential

    """
    Potential for Frohlich Polaron
    """

    α :: Float64
    ω :: Float64
    ħ :: Float64
    κ :: Float64

    function FrohlichInteractionPotential(α::Float64, ω::Float64, ħ::Float64, κ :: Float64)
        new(α, ω, ħ, κ)
    end
end


struct MexicanHatPotential <: OneBodyPotential

    """
    Mexican Hat potential for a single body.
    """

    ω :: Float64
    function MexicanHatPotential(ω::Float64)
        new(ω)
    end
end


struct CoulombPotential <: TwoBodyPotential

    """
    Coulomb interaction between two bodies.
    """

    κ :: Float64
    function CoulombPotential(κ::Float64)
        new(κ)
    end
end

#types of estimators
abstract type Estimator end

#Energy estimators

struct SimpleEstimator <: Estimator end #Estimator using basic sum of kinetic and potential total_energy

struct SimpleVirialEstimator <: Estimator end #Estimator using virial theorem for potential term

struct ThermodynamicEstimator <: Estimator end #Estimator using thermodynamic theory

struct VirialEstimator <: Estimator end #Estimator derived using virial theorem