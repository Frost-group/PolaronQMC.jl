# types.jl

using CircularArrays

"""
Regime for the simuation to run in
"""

abstract type Regime end


struct Simple_Regime <: Regime  #simple form of calculating action
    function Simple_Regime()
        new()
    end
end

struct Primitive_Regime <: Regime #Calculating using the primitive approximation as per Ceperly paper
    function Primitive_Regime()
        new()
    end
end

"""
Adjuster type to container information about shift width and allow for its auto adjustment

"""

abstract type Adjuster end


#Adjuster for the Single! move algorithm
mutable struct Single_Adjuster <: Adjuster
    adjust_counter :: Int
    shift_width :: Float64
    adjust_unit :: Float64 #how much shift width is adjusted by each time
    function Single_Adjuster(λ::Float64, τ::Float64)
        shift_width = sqrt(4 * λ * τ) * 0.5
        adjust_unit = shift_width
        new(0,shift_width, adjust_unit)
    end
end


#Adjuster for the Displace! move algorithm
mutable struct Displace_Adjuster <: Adjuster
    adjust_counter :: Int
    shift_width :: Float64 
    adjust_unit :: Float64 #how much shift width is adjusted by each time
    function Displace_Adjuster(λ::Float64, τ::Float64)
        shift_width = sqrt(4 * λ * τ)
        adjust_unit = shift_width
        new(0,shift_width, adjust_unit)
    end
end




#Adjuster for the Bisect! move alogrithm
mutable struct Bisect_Adjuster <: Adjuster
    segment_length :: Int
    max_level :: Int
    function Bisect_Adjuster(n_beads::Int)
        if n_beads < 100 #adjusting number of beads adjusted by Bisect dependent on total number of beads
            segment_length = 9
        elseif n_beads > 100 && n_beads < 3000
            segment_length = 17
        else
            segment_length = 33
        end

        max_level = Int(floor(log(segment_length)/log(2))) 
        new(segment_length, max_level)
    end
end


"""
Generic path mutable type 
"""
mutable struct Path
	n_beads :: Int64
	n_particles :: Int64
    n_dimensions :: Int64

	beads :: CircularArray{Float64, 3}
    adjusters :: Dict

	τ :: Float64
    m :: Float64
	λ :: Float64



	function Path(n_beads::Int64, n_particles::Int64, n_dimensions::Int64, τ::Float64; m = 1.0, λ = 0.5, start_range = 1.0)
        beads = CircularArray(rand(n_beads, n_particles, n_dimensions) .* (rand([-1,1] * start_range, n_beads, n_particles, n_dimensions)))

        #creating adjusters
        adjusters = Dict()
        adjusters["Single!"] = Single_Adjuster(λ,τ)
        adjusters["Displace!"] = Displace_Adjuster(λ,τ)
        adjusters["Bisect!"] = Bisect_Adjuster(n_beads)



		new(n_beads, n_particles, n_dimensions, beads, adjusters, τ, m, λ)
	end
end

"""
Defined Types for different potentials.
"""

# Most abstract Potential type.
abstract type Potential end

# A potential that is independent on bodies.
abstract type NoBodyPotential <: Potential end

# For potentials that depend on one body.
abstract type OneBodyPotential <: Potential end

# For potentials that depend on two bodies.
abstract type TwoBodyPotential <: Potential end

# A constant potential.
struct ConstantPotential <: NoBodyPotential 
    V :: Float64
    function ConstantPotential(V::Float64)
        new(V)
    end
end

# A Harmonic potential for a single body.
struct HarmonicPotential <: OneBodyPotential
    ω :: Float64
    function HarmonicPotential(ω::Float64)
        new(ω)
    end
end

struct FrohlichPotential <: OneBodyPotential
    α :: Float64
    ω :: Float64
    ħ :: Float64
    function FrohlichPotential(α::Float64, ω::Float64, ħ::Float64)
        new(α, ω, ħ)
    end
end

# Mexican Hat potential for a single body.
struct MexicanHatPotential <: OneBodyPotential
    ω :: Float64
    function MexicanHatPotential(ω::Float64)
        new(ω)
    end
end

# Coulomb interaction between two bodies.
struct CoulombPotential <: TwoBodyPotential
    κ :: Float64
    function CoulombPotential(κ::Float64)
        new(κ)
    end
end







