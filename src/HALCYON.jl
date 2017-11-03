# HALCYON.jl
# Pigs might just fly.
# Path Integral Density Inference Theory

module HALCYON

# Namespace is currently very dirty and everything will need to get (re)sorted
# as it becomes clear where things should lie.
include("system.jl")
include("potentials.jl") #1D potentials for testing
include("Schrodinger.jl") # discretised TISE solver
include("PIMC.jl") # Path integrals for the win

end # module
