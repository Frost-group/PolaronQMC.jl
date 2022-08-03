begin
using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility
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
    n_beads = 80
    τ = 1.0 / (T * n_beads)
    n_particles = 1
    start_range = 1.0
    path = Path(n_beads, n_particles, τ, m = m, λ = λ, start_range = start_range)


    #for comparison energy
    ħ = 1.0
    β = 1/T
    T_scale_factor = 1/7.61

    #for pimc
    n_steps = 1000
    equilibrium_skip = 0.1*n_steps
    #equilibrium_skip = 0
    observables_skip = 0.01*n_steps
    #observables_skip = 100
    movers = [[Single!, Displace!],[1.0, 0.2]]
    #adjusters = [Bisect_Adjuster(path),Single_Adjuster(path)]

    observables = [Energy]
    
    
    #potential type
        potential = FrohlichPotential(α,ω,ħ,β)
        #potential = HarmonicPotential(ω)
    
    #estimator type
    estimators = [Virial_Estimator(n_beads)]
        #estimator = Simple_Estimator()
        #estimators = [Thermodynamic_Estimator()]
        #estimator = Virial_Estimator(n_beads)

    #regime type
        regime = Primitive_Regime()



#running sim


pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true)
acceptance_ratio = pimc[1]
output_observables = pimc[2]

energy = output_observables["Energy"][string(Symbol(estimators[1]))]
#analytic_energy = analytic_energy_harmonic(potential,β,ħ)
comparison_polaron = make_polaron([α], [T*T_scale_factor], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
variances = jackknife(energy)


#post analysis

println("acceptance ratio = ", acceptance_ratio)
println("Mean energy = ", mean(energy))
#println("analytic energy = ", analytic_energy)
println("comparison_energy = ", -comparison_polaron.F)
println("jackknife errors = ", sqrt(variances[2]))



#Plots
worldline = plot(path.beads[:,1], 1:path.n_beads)
energyplot = plot(energy)

plot(worldline, energyplot, layout = (2,1), legend = false)








end

