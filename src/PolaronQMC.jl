# PolaronQMC.jl

module PolaronQMC

using Statistics
using LinearAlgebra

# Export types
export Path, Potential, ConstantPotential, OneBodyPotential, TwoBodyPotential, HarmonicPotential, MexicanHatPotential, CoulombPotential

# Export potentials
export one_body_potential, two_body_potential

# Export observables
export potential_energy, kinetic_energy, Energy, Correlation

# Export PIMC move algorithms
export Single!, Displace!, Bisect!, Stage!

# Export PIMC algorithm
export potential_action, kinetic_action, primitive_action, PIMC, draw_beads_3d, animate_PIMC

include("types.jl") # Types. 
include("potentials.jl") # Potentials for QMC algorithms.
include("moves.jl") # PIMC moves.
include("observables.jl") # Observables to sample.
include("PIMC.jl") # Path integrals Monte Carlo algorithm.
include("Schrodinger.jl") # discretised TISE solver

end # module

