# PolaronQMC.jl

module PolaronQMC

using Revise
using Statistics
using LinearAlgebra
using OnlineStats
using Plots

# Export types
export Path, Potential, ConstantPotential, OneBodyPotential, TwoBodyPotential, HarmonicPotential, MexicanHatPotential, CoulombPotential, Regime, Simple_Regime, Primitive_Regime

# Export potentials
export one_body_potential, two_body_potential

# Export estimators
export potential_energy, kinetic_energy, Energy, Correlation, Estimator, Thermodynamic_Estimator, Virial_Estimator, Simple_Estimator, Simple_Virial_Estimator

# Export PIMC move algorithms
export Single!, Displace!, Bisect!, Stage!

# Export actions
export potential_action, kinetic_action, primitive_action

# Export simulation run code (PIMC)
export PIMC

# Export values for comparison
export analytic_energy_harmonic

# Export error analysis
export jackknife




# Export PIMC algorithm
#export PIMC, draw_beads_3d, animate_PIMC

include("types.jl") # Types. 
include("comparison.jl") #Comparison with actual/experimental values
include("actions.jl") # Kinetic and potential actions.
include("potentials.jl") # Potentials for QMC algorithms.
include("moves.jl") # PIMC moves.
include("estimators.jl") # Estimators to sample.
include("PIMC.jl") # Path integrals Monte Carlo algorithm.
include("errors.jl") # Methods for analysing errors.
#include("Schrodinger.jl") # discretised TISE solver

end # module

