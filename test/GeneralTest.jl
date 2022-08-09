begin
using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility
end


#initialising variables
begin
    #for Potential
    ω = 1.0
    α = 4.0

    #for path
    T = 1.0
    λ = 0.5
    m = ω
    n_beads = 100
    τ = 1.0 / (T * n_beads)
    n_particles = 1
    n_dimensions = 3
    start_range = 1.0
    path = Path(n_beads, n_particles, n_dimensions, τ, m = m, λ = λ, start_range = start_range)


    


    #for comparison energy
    ħ = 1.0
    β = 1/T

    #for pimc
    n_steps = 30000
    #equilibrium_skip = 0.2*n_steps
    equilibrium_skip = 0
    observables_skip = 0.01*n_steps
    #observables_skip = 100
    #movers = [[Bisect!],[1.0]]
    movers = [[Single!],[1.0]]
    #movers = [[Single!, Displace!],[1.0, 0.3]]
    #movers = [[Bisect!, Displace!, Single!],[1.0,0.2,1.0]]

    observables = [Energy, Position]
    
    
    #potential type
        potential = FrohlichPotential(α,ω,ħ,β)
        #potential = HarmonicPotential(ω)
    
    #estimator type
    estimators = [Virial_Estimator(n_beads)]
        #estimator = [Simple_Estimator()]
        #estimators = [Thermodynamic_Estimator()]
        
    #regime type
        regime = Primitive_Regime()



#running sim


pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true)
acceptance_ratio = pimc[1]
output_observables = pimc[2]

energy = output_observables["Energy"][string(Symbol(estimators[1]))]
position = output_observables["Position"][string(Symbol(estimators[1]))]

#=
# comparison energy
    if typeof(potential) == HarmonicPotential
        comparison_energy = analytic_energy_harmonic(potential,β,ħ)
    elseif typeof(potential) == FrohlichPotential
        comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
        comparison_energy = -comparison_polaron.F
    end
=#

variances = jackknife(energy)

#post analysis


println("acceptance ratio = ", acceptance_ratio)
println("Mean energy = ", mean(energy))
#println("comparison_energy = ", comparison_energy)
#println("Mean Position = ", mean(position))
println("jackknife errors = ", sqrt(variances[2]))



#Plots
#posplot = plot(position, ylabel="Mean position")
energyplot = plot(energy, ylabel="Mean energy", xlabel="n_steps")
posplot = histogram(position)
plot(posplot, energyplot, layout = (2,1), legend = false)




end

#More plotting

begin
    plot([path.beads[1]],[path.beads[2]],[path.beads[3]])
end



