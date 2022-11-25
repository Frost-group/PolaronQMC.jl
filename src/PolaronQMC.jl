# PolaronQMC.jl

module PolaronQMC

using Revise
using Statistics
using LinearAlgebra
using OnlineStats
using Plots
#temporary

# Export types
export Path, Potential, ConstantPotential, OneBodyPotential, TwoBodyPotential, HarmonicPotential, MexicanHatPotential, CoulombPotential, FrohlichPotential, Regime, Simple_Regime, Primitive_Regime

# Export potentials
export one_body_potential, one_body_potential_convert, two_body_potential

# Export observables
export potential_energy, kinetic_energy, Energy, Correlation, Position, Mean_Position

# Export estimators
export Estimator, Thermodynamic_Estimator, Virial_Estimator, Simple_Estimator, Simple_Virial_Estimator, Virial_EstimatorX

# Export adjusters for automatically changing shift width
export Adjuster, Single_Adjuster, Displace_Adjuster, Bisect_Adjuster

# Export function used in the optimisation of the running of the simulation
export update_shift_width!, thermalised_start!

# Export PIMC move algorithms
export Single!, Displace!, Bisect!, Stage!

# Export actions
export potential_action, kinetic_action

# Export simulation run code (PIMC)
export PIMC, PIMCX

# Export methods for comparison
export analytic_energy_harmonic, selective_mean

# Export error analysis
export jackknife

# Export visualisation (temp)
export draw_beads_3d, draw_beads_2d, animate_PIMC

# Export quick runs
export quickrun_frohlich


# Export PIMC algorithm
#export PIMC, draw_beads_3d, animate_PIMC

include("types.jl") # Types. 
include("comparison.jl") #Comparison with actual/experimental values
include("actions.jl") # Kinetic and potential actions.
include("potentials.jl") # Potentials for QMC algorithms.
include("moves.jl") # PIMC moves.
include("optimisation.jl") # Adjusters to shift width used in moves.
include("estimators.jl") # Estimators to sample.
include("PIMC.jl") # Path integrals Monte Carlo algorithm.
include("errors.jl") # Methods for analysing errors.
include("quickruns.jl") #Functions for obtaining results using default parameters.


#temporary
include("PolaronQMCVisualisation.jl")
#include("Schrodinger.jl") # discretised TISE solver

end # module

