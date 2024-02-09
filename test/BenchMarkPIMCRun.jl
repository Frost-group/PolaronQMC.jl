begin
    using Revise
    using PolaronQMC
    using Statistics
    using Plots
    using PolaronMobility
    using LaTeXStrings
    using BenchmarkTools
    include("GeneralPIMC.jl")
end

@btime begin
    n_steps = 200000
    mean_energy,
    jacknife_errors,
    comparison_energy,
    energies,
    positions,
    correlations,
    acceptance_rates,
    adjuster_values,
    equilibrium_skip,
    observables_skip,
    version,
    potential,
    n_beads,
    data = generalPIMC(
        0.1, #Temperature
        1.0, # mass
        1.0, # ω (has to be float)
        4.0, # α (has to be float)
        1, # no of particles
        3, # number of n_dimensions
        "Primitive", # regime type
        true, # fixing beads or not
        0.03, # fixed_τ
        80, # n_beads (if τ is not fixed)
        n_steps, # No. of steps
        50000, # number of thermalisation steps
        "Single", # movers
        "Harmonic", # potential type
        "Virial", # estimators
        false, # Not quick steps
        false, # Not threading
        1.0, # Start Range
        1,
        1,
        0.01, # observable skips
        0.2, # equilibrium skips
        1,
        true,
        1,
    )
end
