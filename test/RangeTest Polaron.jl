#testing for the agreeance of the simulation over a range
begin
using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility
end


begin

#other variables
λ = 0.5
n_particles = 1
start_range = 1.0
ħ = 1.0
n_beads = 100

#for pimc
n_steps = 100000
equilibrium_skip = 0.2*n_steps
observables_skip = 0.01*n_steps
movers = [[Single!],[1.0]]
observables = [Energy]
regime = Primitive_Regime()




end

begin
#temperature --------------------------------------------

    #kept constant
    α = 1.0
    ω = 0.5
    m = ω


    #Estimators and potentials
        estimators = [Thermodynamic_Estimator()]
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
    #scatter!(temp_range,observables_range_T[string(Symbol(estimators[2]))], yerr = errors_range_T[string(Symbol(estimators[2]))], label=string(Symbol(estimators[2])))
    
end



#Testing α range ------------------------------------------------
begin
    #kept constant
    T_scale_factor = 1/7.61

    T = 1.0
    τ = 1.0 / (T * n_beads)
    ħ = 1.0
    β = 1/T
    ω = 1.0
    m = ω

    #Estimators and potentials
        estimators = [Thermodynamic_Estimator()]

    #changing
    alpha_range = 0.0:20.0
    comparison_polaron = make_polaron(alpha_range, [T*T_scale_factor], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)


    #output

    observables_range_L = Dict()
    for estimator in estimators
        observables_range_L[string(Symbol(estimator))] = []
    end

    errors_range_L = Dict()
    for estimator in estimators
        errors_range_L[string(Symbol(estimator))] = []
    end

    comparison_range_L = -1*comparison_polaron.F

    #Starting simulation
    for L in alpha_range
        α = L
        potential = FrohlichPotential(α,ω,ħ,β)

        path = Path(n_beads, n_particles, τ, m = m, λ = λ, start_range = start_range)
        pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime)

        output_observables = pimc[2]
        energy = output_observables["Energy"]

        for estimator in estimators
            estimator_energy = energy[string(Symbol(estimator))]
            
            #energy
            append!(observables_range_L[string(Symbol(estimator))],mean(estimator_energy))

            #errors
            error = sqrt(jackknife(estimator_energy)[2])
            append!(errors_range_L[string(Symbol(estimator))], error)

        end
    end



    #Plotting --------------
    plot(alpha_range,comparison_range_L, label="Comparison Energy",xlabel="α",ylabel="Energy",linestyle=:dash,linecolor=:red, linewidth = 1.5)

    #Auto plotting all estimators 
    #=
    for estimator in estimators
        scatter!(temp_range,observables_range_T[string(Symbol(estimator))], yerr = errors_range_T[string(Symbol(estimator))], label=string(Symbol(estimator)))
    end
    =#


    scatter!(alpha_range,observables_range_L[string(Symbol(estimators[1]))], yerr = errors_range_L[string(Symbol(estimators[1]))], label=string(Symbol(estimators[1])))
    #scatter!(alpha_range,observables_range_L[string(Symbol(estimators[2]))], yerr = errors_range_L[string(Symbol(estimators[2]))], label=string(Symbol(estimators[2])))

end


end #final end