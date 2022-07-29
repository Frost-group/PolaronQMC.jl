#testing for the agreeance of the simulation over a range
begin
using Revise
using PolaronQMC
using Statistics
using Plots
end


begin

#other variables
λ = 0.5
n_particles = 1
start_range = 1.0
ħ = 1.0
n_beads = 100

#for pimc
n_steps = 200000
equilibrium_skip = 0.1*n_steps
observables_skip = 0.01*n_steps
movers = [[Single!,Displace!],[1.0,0.1]]
observables = [Energy]
regime = Primitive_Regime()





#temperature --------------------------------------------

#kept constant
ω = 0.5
m = ω

#Estimators and potentials

    #estimator = Simple_Estimator()
    #estimator = Thermodynamic_Estimator()
    estimator = Virial_Estimator(100)

    potential = HarmonicPotential(ω)


#changing
temp_range = 1:10
observables_range = []
errors_range = []
comparison_range = []


for T in temp_range
    τ = 1.0 / (T * n_beads)
    β = 1/T

    path = Path(n_beads, n_particles, τ, m = m, λ = λ, start_range = start_range)
    pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimator, potential, regime)

    output_observables = pimc[2]
    energy = output_observables["Energy"]
    analytic_energy = analytic_energy_harmonic(potential,β,ħ)
    error = sqrt(jackknife(energy)[2])



    append!(observables_range,mean(energy))
    append!(errors_range, error)
    append!(comparison_range, analytic_energy)
    
end

plot(temp_range,comparison_range, label="Analytic Energy",xlabel="Temperature",ylabel="Energy",linestyle=:dash,linecolor=:red, linewidth = 1.5)
scatter!(temp_range,observables_range, yerr = errors_range, label="Simulated Energy")

end

