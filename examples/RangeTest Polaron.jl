#testing for the agreeance of the simulation over a range
begin
using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility
end


#Thermalised start + Alpha test
begin

    st_n_particles = 1
    st_start_range = 1.0

    st_T = 3.0
    st_n_steps = 80000
    st_n_dimensions = 3
    #ability to change beads based of temperature instead.
    fixed_τ = 0.002
    adjusted_beads = Int(floor(1.0 / (st_T * fixed_τ)))
    println("Number of beads = ",adjusted_beads)

    st_n_beads = 50
    st_τ = 1.0 / (st_T * st_n_beads)
    
    #for starting potential
    st_ω = 1.0
    st_α = 1.0
    st_ħ = 1.0



    #path = Path(st_n_beads, st_n_particles, st_n_dimensions, st_τ)
    path = Path(adjusted_beads, st_n_particles, st_n_dimensions, fixed_τ)
    st_movers = [[Single!],[1.0]]

    st_potential = FrohlichPotential(st_α,st_ω,st_ħ)
    #thermalised_start!(path, st_potential, n_steps = st_n_steps)


    #Testing α range ------------------------------------------------

    #kept constant
    n_particles = 1
    start_range = 1.0
    ħ = 1.0
    n_dimensions = 3


    n_beads = 25
    T = st_T
    τ = 1.0 / (T * n_beads)



    ħ = 1.0
    β = 1/T
    ω = st_ω
    m = ω

    #Estimators, movers and other components of PIMC
    estimators = [Virial_Estimator()]
    #estimators = [Thermodynamic_Estimator()]
    movers = [[Bisect!],[1.0]]
    observables = [Energy]
    regime = Primitive_Regime()
        

    #changing
    alpha_range = 1.0:1.0:4.0
    comparison_polaron = make_polaron(alpha_range, [T], [0.0]; ω = 1.0, rtol = 1e-4, verbose = true, threads = true)


    #output

    observables_range_L = Dict()
    for estimator in estimators
        observables_range_L[string(Symbol(estimator))] = []
    end

    errors_range_L = Dict()
    for estimator in estimators
        errors_range_L[string(Symbol(estimator))] = []
    end

    comparison_range_L = comparison_polaron.F

    #Starting simulation
    for L in alpha_range
        
        α = L
        n_steps = Int(10000 + 10000*L)
        #n_steps = 20000
        potential = FrohlichPotential(α,ω,ħ)
        equilibrium_skip = 0.7*n_steps
        observables_skip = 0.03*n_steps

        path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)
        pimc = PIMCX(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true)

        output_observables = pimc[2]
        energy = output_observables["Energy"]

        for estimator in estimators
            estimator_energy = energy[string(Symbol(estimator))]
            
            #energy
            append!(observables_range_L[string(Symbol(estimator))], mean(estimator_energy))

            #errors
            error = sqrt(jackknife(estimator_energy)[2])
            append!(errors_range_L[string(Symbol(estimator))], error)

        end
    end

    println("term one =", 3*T/2)

    #Plotting --------------
    plot(alpha_range,comparison_range_L, label="Comparison Energy",xlabel="α",ylabel="Energy",linestyle=:dash,linecolor=:red, linewidth = 1.5, legend=:bottomleft)

    #Auto plotting all estimators 
    #=
    for estimator in estimators
        scatter!(temp_range,observables_range_T[string(Symbol(estimator))], yerr = errors_range_T[string(Symbol(estimator))], label=string(Symbol(estimator)))
    end
    =#


    scatter!(alpha_range,observables_range_L[string(Symbol(estimators[1]))], yerr = errors_range_L[string(Symbol(estimators[1]))], label=string(Symbol(estimators[1])))
    #scatter!(alpha_range,observables_range_L[string(Symbol(estimators[2]))], yerr = errors_range_L[string(Symbol(estimators[2]))], label=string(Symbol(estimators[2])))

end



begin
addon = 0.5 * coth(β)

    plot(alpha_range,comparison_range_L, label="Comparison Energy",xlabel="α",ylabel="Energy",linestyle=:dash,linecolor=:red, linewidth = 1.5, legend=:bottomleft )

    scatter!(alpha_range,observables_range_L[string(Symbol(estimators[1]))] .+ addon, yerr = errors_range_L[string(Symbol(estimators[1]))], label=string(Symbol(estimators[1])))
    #scatter!(alpha_range,observables_range_L[string(Symbol(estimators[1]))] .+ addition, yerr = errors_range_L[string(Symbol(estimators[2]))], label=string(Symbol(estimators[2])))
end



# beads test  ---------

begin

    #kept constant
    n_particles = 1
    start_range = 1.0
    ħ = 1.0
    n_dimensions = 3


    T = 10.0
    ħ = 1.0
    ω = 1.0
    m = 1.0
    α = 2.0

    #Estimators, movers and other components of PIMC
    estimators = [Virial_EstimatorX()]
    #estimators = [Thermodynamic_Estimator()]
    movers = [[Single!], Displace!,[1.0, 0.2]]
    observables = [Energy]
    regime = Primitive_Regime()
        

    #changing
    beads_step = 50
    beads_range = 1:10


    #output

    observables_range_L = Dict()
    for estimator in estimators
        observables_range_L[string(Symbol(estimator))] = []
    end



    #Starting simulation
    for L in beads_range

        
        n_beads = L*beads_step
        n_steps = 30000*L
        τ = 1.0 / (T * n_beads)

        potential = FrohlichPotential(α,ω,ħ)
        equilibrium_skip = 0.5*n_steps
        observables_skip = 0.05*n_steps

        path = Path(n_beads, n_particles, n_dimensions, τ)
        pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true)

        output_observables = pimc[2]
        energy = output_observables["Energy"]

        for estimator in estimators
            estimator_energy = energy[string(Symbol(estimator))]
            
            #energy
            append!(observables_range_L[string(Symbol(estimator))], mean(estimator_energy))

        end
    end



    #Plotting --------------

    plot(beads_range*beads_step,observables_range_L[string(Symbol(estimators[1]))], label=string(Symbol(estimators[1])), xlabel="Number of beads", ylabel="Energy")
    scatter!(beads_range*beads_step,observables_range_L[string(Symbol(estimators[1]))])
    #scatter!(alpha_range,observables_range_L[string(Symbol(estimators[2]))], yerr = errors_range_L[string(Symbol(estimators[2]))], label=string(Symbol(estimators[2])))

end



#Temperature ---------------------------------------------
begin
    
    #kept constant
    n_particles = 1
    start_range = 1.0
    n_dimensions = 3

    fixed_τ = 0.01



    α = 1.0
    ħ = 1.0
    ω = 1.0
    m = ω

    #Estimators, movers and other components of PIMC
    estimators = [Virial_EstimatorX()]
    #estimators = [Thermodynamic_Estimator()]
    movers = [[Bisect!],[1.0]]
    observables = [Energy]
    regime = Primitive_Regime()
        

    #changing
    T_range = 1.0:1.0:4.0
    comparison_polaron = make_polaron([α], T_range, [0.0]; ω = 1.0, rtol = 1e-4, verbose = true, threads = true)


    #output

    observables_range_T = Dict()
    for estimator in estimators
        observables_range_T[string(Symbol(estimator))] = []
    end

    errors_range_T = Dict()
    for estimator in estimators
        errors_range_T[string(Symbol(estimator))] = []
    end

    comparison_range_T = reshape(comparison_polaron.F,length(T_range))

    #Starting simulation
    for T in T_range
        
        n_steps = Int(7000 + (2000 * T))
        adjusted_beads = Int(floor(1.0 / (T * fixed_τ)))
        

        #n_steps = 20000
        potential = FrohlichPotential(α,ω,ħ)
        equilibrium_skip = 0.7*n_steps
        observables_skip = 0.03*n_steps

        path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)
        pimc = PIMCX(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true)

        output_observables = pimc[2]
        energy = output_observables["Energy"]

        for estimator in estimators
            estimator_energy = energy[string(Symbol(estimator))]
            
            #energy
            append!(observables_range_T[string(Symbol(estimator))], mean(estimator_energy))

            #errors
            error = sqrt(jackknife(estimator_energy)[2])
            append!(errors_range_T[string(Symbol(estimator))], error)

        end
    end


    #Plotting --------------
    plot(T_range,comparison_range_T, label="Comparison Energy",xlabel="T",ylabel="Energy",linestyle=:dash,linecolor=:red, linewidth = 1.5, legend=:bottomleft)

    #Auto plotting all estimators 
    #=
    for estimator in estimators
        scatter!(temp_range,observables_range_T[string(Symbol(estimator))], yerr = errors_range_T[string(Symbol(estimator))], label=string(Symbol(estimator)))
    end
    =#


    scatter!(T_range,observables_range_T[string(Symbol(estimators[1]))], yerr = errors_range_T[string(Symbol(estimators[1]))], label=string(Symbol(estimators[1])))
    #scatter!(alpha_range,observables_range_L[string(Symbol(estimators[2]))], yerr = errors_range_L[string(Symbol(estimators[2]))], label=string(Symbol(estimators[2])))

end
 
