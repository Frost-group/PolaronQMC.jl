begin
using BenchmarkTools
using Revise
using PolaronQMC
end


begin
println("Using ", Threads.nthreads(), " threads")
#initialising
    α = 1.0
    ω = 1.0
    ħ = 1.0
    T = 1.0
    n_beads = 100
    τ = 1.0 / (T * n_beads)
    n_particles = 1
    n_dimensions = 3

    path = Path(n_beads, n_particles, n_dimensions, τ)

#variables --------------------------------------
    potential = FrohlichPotential(α,ω,ħ)

    #timing

    println("Time take for Virial_Estimator:")
    thermalised_start!(path, potential, n_steps = 500)
    energy_bm = @btime Energy(path, potential, Virial_Estimator())

end