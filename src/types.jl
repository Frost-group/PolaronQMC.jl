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
Generic path mutable type 
"""

#One dimensional version


mutable struct Path
	n_beads :: Int64
	n_particles :: Int64

	beads :: CircularArray{Float64, 2}

	τ :: Float64
    m :: Float64
	λ :: Float64


	function Path(n_beads::Int64, n_particles::Int64, τ::Float64; m = 1.0, λ = 0.5, start_range = 1.0)
        beads = CircularArray(rand(n_beads, n_particles).*rand([-1,1],n_beads)*start_range)
		new(n_beads, n_particles, beads, τ, m, λ)
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
    β :: Float64
    function FrohlichPotential(α::Float64, ω::Float64, ħ::Float64, β::Float64)
        new(α, ω, ħ, β)
    end
end

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




"""
Adjuster type to container information about shift width and allow for its auto adjustment

"""


abstract type Adjuster end

mutable struct Single_Adjuster <: Adjuster
    adjust_counter :: Int
    shift_width :: Float64
    adjust_unit :: Float64 #how much shift width is adjusted by each time
    function Single_Adjuster(path::Path)
        shift_width = sqrt(4 * path.λ * path.τ)
        adjust_unit = shift_width
        new(0,shift_width, adjust_unit)
    end
end

mutable struct Displace_Adjuster <: Adjuster
    adjust_counter :: Int
    shift_width :: Float64 
    adjust_unit :: Float64 #how much shift width is adjusted by each time
    function Displace_Adjuster(path::Path)
        shift_width = sqrt(4 * path.λ * path.τ)*100
        adjust_unit = 0.5*shift_width
        new(0,shift_width, adjust_unit)
    end
end

mutable struct Bisect_Adjuster <: Adjuster
    adjust_counter :: Int
    shift_width :: Float64
    #shift_width_array :: Dict
    adjust_unit :: Float64 #how much shift width is adjusted by each time


    function Bisect_Adjuster(path::Path)
        shift_width = path.λ * path.τ
        adjust_counter = 0
        #shift_width_array = Dict()
        #=
        max_level = Int(floor(log(rand(1:path.n_beads)) / log(2)))

        for level in 1:max_level-1
            shift_width[string(level)] = sqrt(2^(level-1) * path.λ * path.τ )
            adjust_counter[string(level)] = 0
        end
        =#
        adjust_unit = 0.5*shift_width
        new(adjust_counter,shift_width, adjust_unit)
    end
end


