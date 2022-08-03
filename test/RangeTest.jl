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
n_steps = 40000
equilibrium_skip = 0.1*n_steps
observables_skip = 0.02*n_steps
movers = [[Single!, Bisect!, Displace!],[1.0, 0.3, 0.2]]
#movers = [[Single!, Displace!],[1.0, 0.2]]
observables = [Energy]

regime = Primitive_Regime()


end

begin
#temperature --------------------------------------------

    #kept constant
    ω = 1.0
    m = ω


    #Estimators and potentials
        estimators = [Virial_Estimator(n_beads)]
        #estimators = [Thermodynamic_Estimator()]
        potential = HarmonicPotential(ω)


    #changing
    temp_range = 1:10


    #output

    observables_range_T = Dict()
    for estimator in estimators
        observables_range_T[string(Symbol(estimator))] = []
    end

    errors_range_T = Dict()
    for estimator in estimators
        errors_range_T[string(Symbol(estimator))] = []
    end

    comparison_range_T = []

    #Starting simulation
    for T in temp_range
        τ = 1.0 / (T * n_beads)
        β = 1/T

        path = Path(n_beads, n_particles, τ, m = m, λ = λ, start_range = start_range)
        pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust = true)

        output_observables = pimc[2]
        analytic_energy = analytic_energy_harmonic(potential,β,ħ)
        energy = output_observables["Energy"]

        for estimator in estimators
            estimator_energy = energy[string(Symbol(estimator))]
            
            #energy
            append!(observables_range_T[string(Symbol(estimator))],mean(estimator_energy))

            #errors
            error = sqrt(jackknife(estimator_energy)[2])
            append!(errors_range_T[string(Symbol(estimator))], error)

        end
        append!(comparison_range_T, analytic_energy)
        
    end

    plot(temp_range,comparison_range_T, label="Analytic Energy",xlabel="Temperature",ylabel="Energy",linestyle=:dash,linecolor=:red, linewidth = 1.5)

    #=
    for estimator in estimators
        scatter!(temp_range,observables_range_T[string(Symbol(estimator))], yerr = errors_range_T[string(Symbol(estimator))], label=string(Symbol(estimator)))
    end
    =#
    
    
    scatter!(temp_range,observables_range_T[string(Symbol(estimators[1]))], yerr = errors_range_T[string(Symbol(estimators[1]))], label=string(Symbol(estimators[1])))
    #scatter!(temp_range,observables_range_T[string(Symbol(estimators[2]))], yerr = errors_range_T[string(Symbol(estimators[2]))], label=string(Symbol(estimators[2])))
    
end

begin
    plot(temp_range,comparison_range_T, label="Analytic Energy",xlabel="Temperature",ylabel="Energy",linestyle=:dash,linecolor=:red, linewidth = 1.5)
    scatter!(temp_range,observables_range_T[string(Symbol(estimators[1]))], yerr = errors_range_T[string(Symbol(estimators[1]))], label=string(Symbol(estimators[1])))
end

#lattice spacing range (ω and m) ------------------------------------------------
begin
    #kept constant
    T = 1.0
    τ = 1.0 / (T * n_beads)
    β = 1/T
    ω = 1.0
    m = ω

    #Estimators and potentials
        estimators = [Thermodynamic_Estimator()]

    #changing
    lattice_range = 1.0:12.0


    #output

    observables_range_L = Dict()
    for estimator in estimators
        observables_range_L[string(Symbol(estimator))] = []
    end

    errors_range_L = Dict()
    for estimator in estimators
        errors_range_L[string(Symbol(estimator))] = []
    end

    comparison_range_L = []

    #Starting simulation
    for L in lattice_range
        α = L
        potential = FrohlichPotential(α,ω)

        path = Path(n_beads, n_particles, τ, m = m, λ = λ, start_range = start_range)
        adjusters = [Single_Adjuster(path)]
        pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjusters)

        output_observables = pimc[2]
        #analytic_energy = analytic_energy_harmonic(potential,β,ħ)
        energy = output_observables["Energy"]

        for estimator in estimators
            estimator_energy = energy[string(Symbol(estimator))]
            
            #energy
            append!(observables_range_L[string(Symbol(estimator))],mean(estimator_energy))

            #errors
            error = sqrt(jackknife(estimator_energy)[2])
            append!(errors_range_L[string(Symbol(estimator))], error)

        end
        #append!(comparison_range_L, analytic_energy)
    end

    #plot(lattice_range,comparison_range_L, label="Analytic Energy",xlabel="m and ω",ylabel="Energy",linestyle=:dash,linecolor=:red, linewidth = 1.5)

    #Auto plotting all estimators 
    #=
    for estimator in estimators
        scatter!(temp_range,observables_range_T[string(Symbol(estimator))], yerr = errors_range_T[string(Symbol(estimator))], label=string(Symbol(estimator)))
    end
    =#
    
    
    plot(lattice_range,observables_range_L[string(Symbol(estimators[1]))], yerr = errors_range_L[string(Symbol(estimators[1]))], label=string(Symbol(estimators[1])))
    #scatter!(lattice_range,observables_range_L[string(Symbol(estimators[2]))], yerr = errors_range_L[string(Symbol(estimators[2]))], label=string(Symbol(estimators[2])))
    


end
end