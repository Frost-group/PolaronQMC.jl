# PolaronQMC.jl

module PolaronQMC

using Revise
using Statistics
using LinearAlgebra
using OnlineStats
using Plots
using JLD
using TimerOutputs
#temporary

# Export types
export Path, DiscretePath, Potential, ConstantPotential, OneBodyPotential, TwoBodyPotential, HarmonicPotential, MexicanHatPotential, CoulombPotential, FrohlichPotential, Regime, SimpleRegime, PrimitiveRegime, BoundRegime, LBRegime, Mover, SingleMover, BisectMover, DisplaceMover, HarmonicInteractionPotential, FrohlichInteractionPotential, HolsteinPotential

# Export adjusters for automatically changing shift width
export Adjuster, SingleAdjuster, DisplaceAdjuster, BisectAdjuster

# Export estimators
export Estimator, ThermodynamicEstimator, VirialEstimator, SimpleEstimator, SimpleVirialEstimator

# Export potentials
export oneBodyPotential, twoBodyPotential

# Export observables
export potentialEnergy, kineticEnergy, energy, correlation, position, HolsteinEnergy

# Export function used in the optimisation of the running of the simulation
export updateAdjuster, copyLastPath!

# Export PIMC move algorithms
export moveBead!

# Export actions
export totalAction, potentialAction, kineticAction, bisectPotentialAction

# Export simulation run code (PIMC)
export PIMC, Holstein_PIMC

# Export methods for comparison
export analyticEnergyHarmonic, selectiveMean

# Export error analysis
export jackknife, autoCorrelation, autoCorrelationTime

# Export visualisation (temp)
export draw_beads_3d, draw_beads_2d, animate_PIMC

# Export quick runs
export quickrun_frohlich

export generalPIMC, MultiModePIMC, general_Holstein_PIMC, SaveJLDData


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
include("quickruns.jl") # Functions for obtaining results using default parameters.
include("GeneralPIMC.jl") # Function for making data


#temporary
include("PolaronQMCVisualisation.jl")

end # module

