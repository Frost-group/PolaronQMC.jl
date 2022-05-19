# types.jl

"""
Generic path mutable type
"""
mutable struct Path
	n_beads :: Int64
	n_particles :: Int64
	n_dimensions :: Int64

	beads :: Array{Float64, 3}

	τ :: Float64
	λ :: Float64

	function Path(n_beads::Int64, n_particles::Int64; n_dimensions::Int64 = 3, τ = 0.05, λ = 0.5)
		beads = rand(n_beads, n_particles, n_dimensions)
		new(n_beads, n_particles, n_dimensions, beads, τ, λ)
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

# Coulomb interaction between two bodies.
struct CoulombPotential <: TwoBodyPotential
    κ :: Float64
    function CoulombPotential(κ::Float64)
        new(κ)
    end
end