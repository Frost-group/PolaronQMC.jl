begin
    using LinearAlgebra
    using Revise
    using PolaronQMC
    using Statistics
    using Plots
    #include("../src/PolaronQMCVisualisation.jl")
    #using .PolaronQMCVisualisation
    using PolaronMobility

end



begin
#initialising variables
    #for Potential
        ω = 1.0
        α = 3.0
        ħ = 1.0
    #for path ---------
        #T = 10.0
        T = 0.1
        m = ω
        n_beads = 500 # So far 1000 has good results
        τ = 1.0 / (T * n_beads)
        n_particles = 1
        n_dimensions = 3
        start_range = 1.0
        β = 1/T
        path = Path(n_beads, n_particles, n_dimensions, τ, m = m)
    # regime
        regime = Primitive_Regime()

#variables more subject to change
    #potential type --------------------------
        potential = FrohlichPotential(α,ω,ħ)
        #potential = HarmonicPotential(ω)

    #for pimc --------------------------------------
        #number of steps
            #n_steps = 40000
            n_steps = 1000 #for alpha = 4

        #skipping between sampling
            #equilibrium_skip = 0.2*n_steps
            equilibrium_skip = 0
            observables_skip = 0.01*n_steps
            observables_skip = 10

        #types of moves
            #movers = [[Single!],[1.0]]
            #movers = Dict("Single!" => [1.0])
            movers = Dict("Bisect!" => [1.0])
            #movers = [[Single!],[1.0]]
            #movers = [[Single!, Displace!],[1.0, 0.3]]

        #observables
            observables = [Energy,Position]
    
        #estimator type
            estimators = [Virial_Estimator()]
            #estimators = [Thermodynamic_Estimator()]
        
#running sim

    #thermalised_start!(path,potential,n_steps = 3000)
    pimc = PIMC(n_steps::Int, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true, visual=true)
    #acceptance_ratio = pimc[1] # older version
    adjuster_stats = pimc[1]
    acceptance_rates = adjuster_stats["Single!"]["Acceptance Rate"]
    shift_widths = adjuster_stats["Single!"]["Shift Width"]
    acceptance_rates = adjuster_stats["Bisect!"]["Acceptance Rate"]
    shift_widths = adjuster_stats["Bisect!"]["Shift Width"]
    output_observables = pimc[2]

    energy = output_observables["Energy"][string(Symbol(estimators[1]))]
    position = output_observables["Position"][string(Symbol(estimators[1]))]

    
    # comparison energy
        if typeof(potential) == HarmonicPotential
            comparison_energy = analytic_energy_harmonic(potential,β,ħ, n_dimensions)
        elseif typeof(potential) == FrohlichPotential
            comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
            #comparison_energy = -comparison_polaron.F
            comparison_energy = comparison_polaron.F
        end
    


    variances = jackknife(energy)

    #post analysis
    last_acceptance_rate = last(acceptance_rates)
    mean_acceptance_rate = mean(acceptance_rates)
    std_acceptance_rate = std(acceptance_rates)


    println("Final Acceptance Rate: ", last_acceptance_rate)
    println("Mean Acceptance Rate: ", mean_acceptance_rate, " +/- ", std_acceptance_rate)
    println("Mean energy = ", mean(energy[Int(0.5*length(energy)):end]))
    #println("comparison_energy = ", dot(comparison_energy, 1) - 3/2/β * log(2*pi*β))
    println("comparison_energy = ", comparison_energy)
    println("jackknife errors = ", sqrt(variances[2]))



    #Plots
    energyplot = plot(energy, ylabel="Mean energy", xlabel="n_steps", title="polaron energy \n α=$(α)")
    posplot = histogram(position[:,1,1])
    acc_plot = plot(acceptance_rates, ylabel="accept rate", xlabel="n_steps", title="acceptance_profile \n α=$(α)")
    plot(posplot, energyplot, acc_plot, layout = (3,1), legend = false)
    #savefig("./figs/Froh/α$(α)vsE.png")
    #plot(energyplot)

    #acc_plot = plot(acceptance_rates, ylabel="accept rate", xlabel="n_steps", title="acceptance_profile \n α=$(α)")
    
end

begin
    energyplot = plot(energy, ylabel="Mean energy", xlabel="n_steps", title="polaron energy \n α=$(α)")
    posplot = histogram(position[:,1,1])
    acc_plot = scatter(acceptance_rates[Int(floor(0.9*length(acceptance_rates))):end], ylabel="accept rate", xlabel="n_steps", title="acceptance_profile \n α=$(α)")
    plot(posplot, energyplot, acc_plot, layout = (3,1), legend = false)
end




#=
begin
# Visualise

anim = animate_PIMC(pimc, n_particles)
gif(anim, "saved_plots/anim_output.gif", fps = 10)

end
=#
