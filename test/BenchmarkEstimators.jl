begin
using BenchmarkTools
using Revise
using PolaronQMC
end


begin
#initialising
    ω = 1.0
    α = 10.0
    ħ = 1.0
    T = 1.0
    m = ω
    n_beads = 100
    τ = 1.0 / (T * n_beads)
    n_particles = 1
    n_dimensions = 3
    start_range = 1.0
    β = 1/T

    path = Path(n_beads, n_particles, n_dimensions, τ, m = m)

    observables = [Energy, Position]
    regime = Primitive_Regime()


#variables --------------------------------------
    n_steps = 2000
    equilibrium_skip = 0.1*n_steps
    observables_skip = 0.01*n_steps
    movers = [[Bisect!],[1.0]]
    estimators = [Virial_Estimator()]
    potential = FrohlichPotential(α,ω,ħ)

    #timing
    pimc = @btime PIMC(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true)
    println("Complete")

end