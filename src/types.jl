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
        adjusters["BisectL!"] = Bisect_Adjuster(n_beads)



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

"""
Cache objects containing constant information used in potential calculation to avoid redundant calculations
"""

abstract type Cache end

"""
    PotentialCache(path::Path, potential::Potential)

A cache of constant information used when evaluating the potential energy.

# Arguments
- `path::Path`: collection of all particle imaginary-time paths.
- `potential::Potential`: the potential controlling the system.
"""
mutable struct PotentialCache <: Cache

    potential_prefactor :: Union{Float64, CircularArray{}}
    distance_matrix :: Union{CircularArray{Float64}, Array}
    old_distance_matrix :: Union{CircularArray{Float64}, Array} #storage of the old distance matrix to revert to in the case of a failed move.

    function PotentialCache(path::Path, potential::FrohlichPotential)

        #generating the g factors used in frohlich potential calculation, combination of constants and adjacency information 
            β = path.n_beads * path.τ
            g_prefactor = -0.5 * potential.α * (potential.ħ * potential.ω)^3/2 * sqrt(2*path.m) * csch(potential.ħ * potential.ω * β / 2)

            function generate_g_factors(bead::Int, path::Path, potential::FrohlichPotential)
                return g_prefactor .* [cosh(potential.ω*β * (abs(bead-other_bead)/path.n_beads - 0.5 * potential.ħ)) for other_bead in 1:path.n_beads if other_bead != bead]
            end

            g_factors_array = CircularArray([generate_g_factors(bead, path, potential) for bead in 1:path.n_beads])

        #generating a matrix containing the distance between beads which will be updated after each successful move, currently only supports 1 particle
            particle = 1
            distance_matrix = CircularArray(reduce(hcat, [generate_distances(bead, particle, path) for bead in 1:path.n_beads])) #array containing the distance between each and every bead
            old_distance_matrix = copy(distance_matrix)

        new(g_factors_array, distance_matrix, old_distance_matrix)
    end

    function PotentialCache(path::Path, potential::HarmonicPotential)
            prefactor_1 = 0.5 * path.m * potential.ω^2



        new(prefactor_1, [], [])

    end
end






