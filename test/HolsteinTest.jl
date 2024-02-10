# HolsteinTest.jl - basic Holstein integration test

pot = "Holstein"
T = 0.05
version = rand(1:10000)
energy, errors, data_set = general_Holstein_PIMC(T, 1.0, 2.0, 3, 10_000_000, version = 1)

