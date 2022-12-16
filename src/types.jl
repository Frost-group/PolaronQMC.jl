# types.jl

using CircularArrays


# Regime for the simuation to run in
abstract type Regime end


struct Simple_Regime <: Regime  # Simple form of calculating action
    function Simple_Regime()
        new()
    end
end


struct Primitive_Regime <: Regime # Calculating using the primitive approximation as per Ceperly paper
    function Primitive_Regime()
        new()
    end
end


mutable struct Path

    """
    Generic path mutable type 
    """

	n_beads :: Int64
	n_particles :: Int64
    n_dimensions :: Int64

	beads :: CircularArray{Float64, 3}
    adjusters :: Dict

	τ :: Float64
    m :: Float64
	λ :: Float64

	function Path(n_beads::Int64, n_particles::Int64, n_dimensions::Int64, τ::Float64; m = 1.0, λ = 0.5, start_range = 1.0)
        beads = CircularArray(rand(n_beads, n_particles, n_dimensions) .* (rand([-1, 1] * start_range, n_beads, n_particles, n_dimensions)))

        # DIctionary of Adjusters
        adjusters = Dict()
        adjusters["Single!"] = Single_Adjuster(λ, τ)
        adjusters["Displace!"] = Displace_Adjuster(λ, τ)
        adjusters["Bisect!"] = Bisect_Adjuster(λ, τ)

		new(n_beads, n_particles, n_dimensions, beads, adjusters, τ, m, λ)
	end
end


# Adjuster type to container information about shift width and allow for its auto adjustment
abstract type Adjuster end


# Adjuster for the Single! move algorithm
mutable struct Single_Adjuster <: Adjuster

    """
    Adjuster for the Single! move algorithm
    """

    attempt_counter :: Int
    success_counter :: Int
    shift_width :: Float64
    acceptance_rate :: Float64
    function Single_Adjuster(λ::Float64, τ::Float64)
        shift_width = sqrt(4 * λ * τ) * 0.5
        new(0, 0, shift_width, 0)
    end
end


#Adjuster for the Displace! move algorithm
mutable struct Displace_Adjuster <: Adjuster

    """
    Adjuster for the Displace! move algorithm
    """

    attempt_counter :: Int
    success_counter :: Int
    shift_width :: Float64
    acceptance_rate :: Float64
    function Displace_Adjuster(λ::Float64, τ::Float64)
        shift_width = 1
        new(0, 0, shift_width, 0)
    end
end


#Adjuster for the Bisect! move algorithm
mutable struct Bisect_Adjuster <: Adjuster

    attempt_counter :: Int
    success_counter :: Int
    shift_width :: Float64
    acceptance_rate :: Float64

    function Bisect_Adjuster(λ::Float64, τ::Float64)
        shift_width = sqrt(τ * λ)
        new(0, 0, shift_width, 0)
    end
end


#=
#Adjuster for the Bisect! move alogrithm
mutable struct Bisect_Adjuster <: Adjuster end

    adjust_counter_array :: Dict
    shift_width_array :: Dict
    max_level :: Int


    function Bisect_Adjuster(λ,τ)
        adjust_counter_array = Dict()
        shift_width_array = Dict()
        
        #max_level = Int(floor(log(rand(1:path.n_beads)) / log(2)))
        max_level = 4

        for level in 0:max_level
            shift_width_array[string(level)] = sqrt(2^(level) * λ * τ )
            adjust_counter_array[string(level)] = 0
        end
        
        new(adjust_counter_array,shift_width_array, max_le vel)
    end
end 
=#

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