begin
using Revise
using PolaronQMC
using Statistics
using Plots
end

begin

#initialising variables

    #for Potential
    ω = 1.0
    α = 10.0

    #for path
    T = 1.0
    λ = 0.5
    m = ω
    n_beads = 1000
    τ = 1.0 / (T * n_beads)
    n_particles = 1
    start_range = 1.0
    path = Path(n_beads, n_particles, τ, m = m, λ = λ, start_range = start_range)


    #for analytic energy
    ħ = 1.0
    β = 1/T

    #for pimc
    n_steps = 100000
    #equilibrium_skip = 0.1*n_steps
    equilibrium_skip = 0
    bservables_skip = 0.01*n_steps
    #observables_skip = 1000
    movers = [[Bisect!],[1.0]]
    adjusters = [Bisect_Adjuster(path)]

    observables = [Energy]
    
    
    #potential type
        #potential = FrohlichPotential(α,ω,ħ,β)
        potential = HarmonicPotential(ω)
    
    #estimator type
    estimators = [Virial_Estimator(500)]
        #estimator = Simple_Estimator()
        #estimator = Thermodynamic_Estimator()
        #estimator = Virial_Estimator(100)s

    #regime type
        regime = Primitive_Regime()



#running sim


pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjusters, adjust=true)
acceptance_ratio = pimc[1]
output_observables = pimc[2]

energy = output_observables["Energy"][string(Symbol(estimators[1]))]
#analytic_energy = analytic_energy_harmonic(potential,β,ħ)
#variances = jackknife(energy)


#post analysis

println("acceptance ratio -", acceptance_ratio)
println("Mean energy = ", mean(energy))
#println("analytic energy = ", analytic_energy)
#println("jackknife errors = ", sqrt(variances[2]))



#Plots
worldline = plot(path.beads[:,1], 1:path.n_beads)
energyplot = plot(energy)

plot(worldline, energyplot, layout = (2,1), legend = false)








end

