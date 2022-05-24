# PolaronQMC.jl

module PolaronQMC

using Statistics
using LinearAlgebra
using OnlineStats
using Plots

# Export types
export Path, Potential, ConstantPotential, OneBodyPotential, TwoBodyPotential, HarmonicPotential, MexicanHatPotential, CoulombPotential

# Export potentials
export one_body_potential, two_body_potential

# Export estimators
export potential_energy, kinetic_energy, Energy, Correlation

# Export PIMC move algorithms
export Single!, Displace!, Bisect!, Stage!

# Export actions
export potential_action, kinetic_action, primitive_action

# Export PIMC algorithm
export PIMC, draw_beads_3d, animate_PIMC

include("types.jl") # Types. 
include("actions.jl") # Kinetic and potential actions.
include("potentials.jl") # Potentials for QMC algorithms.
include("moves.jl") # PIMC moves.
include("estimators.jl") # Estimators to sample.
include("PIMC.jl") # Path integrals Monte Carlo algorithm.
include("Schrodinger.jl") # discretised TISE solver

end # module

