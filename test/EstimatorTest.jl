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

    #for path
    T = 1.0
    λ = 0.5
    m = 1.0
    n_beads = 100
    τ = 1.0 / (T * n_beads)
    n_particles = 1
    start_range = 1


    #for analytic energy
    ħ = 1.0
    β = 1/T

    #for pimc
    n_steps = 10000
    equilibrium_skip = 0.1*n_steps
    observables_skip = 0.01*n_steps
    movers = [Single!]
    observables = [Energy]
    
    #potential type
        potential = HarmonicPotential(ω)
    
    #estimator type
        #estimator = Simple_Estimator()
        estimator = Thermodynamic_Estimator()
        #estimator = Virial_Estimator(n_beads)

    #regime type
        regime = Primitive_Regime()



#running sim

path = Path(n_beads, n_particles, τ, m = m, λ = λ, start_range = start_range)
pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimator, potential, regime)
acceptance_ratio = pimc[1]
output_observables = pimc[2]

energy = output_observables["Energy"]
analytic_energy = analytic_energy_harmonic(potential,β,ħ)


#post analysis

println("acceptance ratio -", acceptance_ratio)
println("Mean energy = ", mean(energy))
println("analytic energy = ", analytic_energy)
println("average position = ", mean(path.beads[:,1]))


#Plots
worldline = plot(path.beads[:,1], 1:path.n_beads)
energyplot = plot(energy)

plot(worldline, energyplot, layout = (2,1), legend = false)








end

#Testing new estimators