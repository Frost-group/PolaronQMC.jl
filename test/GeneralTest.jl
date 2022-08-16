begin
using Revise
using PolaronQMC
using Statistics
using Plots
include("../src/PolaronQMCVisualisation.jl")
using .PolaronQMCVisualisation
#using PolaronMobility

end



begin
#initialising variables
    #for Potential
        ω = 1.0
        α = 1.0
        ħ = 1.0
    #for path ---------
        T = 0.05
        m = 1.0
        
        fixed_τ = 0.02
        adjusted_beads = Int(floor(1/(fixed_τ*T)))
        n_beads = 100
        τ = 1.0 / (T * n_beads)
        n_particles = 1
        n_dimensions = 3
        start_range = 1.0
        β = 1/T
        #path = Path(n_beads, n_particles, n_dimensions, τ, m = m)
        path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)
    # regime
        regime = Primitive_Regime()

#variables more subject to change
    #potential type --------------------------
        #potential = FrohlichPotential(α,ω,ħ)
        #potential = HarmonicPotential(ω)
        potential = MexicanHatPotential(80000.0)
        #potential = ConstantPotential(10.0)

    #for pimc --------------------------------------
        #number of steps
            n_steps = 40000

        #skipping between sampling
            equilibrium_skip = 0.1*n_steps
            #equilibrium_skip = 0
            #observables_skip = 0.1*
            observables_skip = 0.02*n_steps

        #types of moves
            #movers = [[Bisect!],[1.0]]
            movers = [[Single!],[1.0]]
            #movers = [[Displace!],[1.0]]

        #observables
            observables = [Position]
    
        #estimator type
            #estimators = [Virial_EstimatorX()]
            estimators = [Thermodynamic_Estimator()]
        
#running sim

    #thermalised_start!(path,potential,n_steps = 3000)
    pimc = PIMC(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true, visual=true)
    acceptance_ratio = pimc[1]
    output_observables = pimc[2]

    #energy = output_observables["Energy"][string(Symbol(estimators[1]))]
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


    #variances = jackknife(energy)

    #post analysis


    println("acceptance ratio = ", acceptance_ratio)
    #println("Mean energy = ", mean(energy))
    #println("comparison_energy = ", comparison_energy)
    #println("jackknife errors = ", sqrt(variances[2]))



    #Plots
    #energyplot = plot(energy, ylabel="Energy", xlabel="n_steps")
    posplot = histogram(position[:,1,1],bins=150)
    #plot(posplot, energyplot, layout = (2,1), legend = false)
    plot(posplot, xlabel="Position", ylabel="Prob Amplitude", legend = false)

end







begin
# Visualise

anim = animate_PIMC(pimc, n_particles)
gif(anim, "saved_plots/anim_output.gif", fps = 8)

end

