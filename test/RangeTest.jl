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


end

begin
#temperature --------------------------------------------

    #kept constant
    ω = 0.5
    m = ω

    #Estimators and potentials
        estimators = [Thermodynamic_Estimator(), Virial_Estimator(1)]
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
        pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime)

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
    scatter!(temp_range,observables_range_T[string(Symbol(estimators[2]))], yerr = errors_range_T[string(Symbol(estimators[2]))], label=string(Symbol(estimators[2])))
    
end



#lattice spacing range (ω and m) ------------------------------------------------
begin
    #kept constant
    T = 1.0
    τ = 1.0 / (T * n_beads)
    β = 1/T


    #Estimators and potentials
        estimators = [Thermodynamic_Estimator(), Virial_Estimator(100)]

    #changing
    lattice_range = 1.0:10.0


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
        ω = L
        m = L
        potential = HarmonicPotential(ω)

        path = Path(n_beads, n_particles, τ, m = m, λ = λ, start_range = start_range)
        pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime)

        output_observables = pimc[2]
        analytic_energy = analytic_energy_harmonic(potential,β,ħ)
        energy = output_observables["Energy"]

        for estimator in estimators
            estimator_energy = energy[string(Symbol(estimator))]
            
            #energy
            append!(observables_range_L[string(Symbol(estimator))],mean(estimator_energy))

            #errors
            error = sqrt(jackknife(estimator_energy)[2])
            append!(errors_range_L[string(Symbol(estimator))], error)

        end
        append!(comparison_range_L, analytic_energy)
        
    end

    plot(lattice_range,comparison_range_L, label="Analytic Energy",xlabel="m and ω",ylabel="Energy",linestyle=:dash,linecolor=:red, linewidth = 1.5)

    #=
    for estimator in estimators
        scatter!(temp_range,observables_range_T[string(Symbol(estimator))], yerr = errors_range_T[string(Symbol(estimator))], label=string(Symbol(estimator)))
    end
    =#
    
    
    scatter!(lattice_range,observables_range_L[string(Symbol(estimators[1]))], yerr = errors_range_L[string(Symbol(estimators[1]))], label=string(Symbol(estimators[1])))
    scatter!(lattice_range,observables_range_L[string(Symbol(estimators[2]))], yerr = errors_range_L[string(Symbol(estimators[2]))], label=string(Symbol(estimators[2])))
    


end
end