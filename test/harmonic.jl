# Harmonic.jl - basic Harmonic PIMC 

pot = "Harmonic"
T = 0.05
energy, errors, data_set = general_Holstein_PIMC(T, 1.0, 2.0, 3, 1_000_000, version = 1)

