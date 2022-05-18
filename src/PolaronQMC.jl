# PolaronQMC.jl

module PolaronQMC

using Statistics
using LinearAlgebra

# Export types
export Path, Potential, ZeroPotential, OneBodyPotential, TwoBodyPotential, HarmonicPotential, CoulombPotential

# Export potentials
export one_body_potential, two_body_potential

# Export observables
export potential_energy, kinetic_energy

# Export PIMC move algorithms
export Single, Displace, Bisection, Staging

# Export PIMC algorithm
export potential_action, kinetic_action, primitive_action, PIMC

include("types.jl") # Types.
include("potentials.jl") # Potentials for QMC algorithms.
include("moves.jl") # PIMC moves.
include("observables.jl") # Observables to sample.
include("PIMC.jl") # Path integrals Monte Carlo algorithm.
include("Schrodinger.jl") # discretised TISE solver

end # module
