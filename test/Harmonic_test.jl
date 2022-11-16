using Revise
using PolaronQMC
using Statistics
using Plots
include("../src/PolaronQMCVisualisation.jl")
using .PolaronQMCVisualisation
using PolaronMobility
using LaTeXStrings


begin
    
    """
    Initialise System Variables
    """

    # Path variables
    T = 1
    T_arr = 10 .^ range(-1, 0, length = 5)
    m = 1.0
    n_particles = 1
    n_dimensions = 1
    start_range = 1.0
    β = 1 / T
    Mean_energy_arr = []
    Comparison_energy_arr = []
    Error_arr = []

    for T in T_arr
        # For fixed τ 
        fixed_τ = 0.01
        adjusted_beads = Int(floor(1/(fixed_τ*T)))

        # For fixed number of beads
        n_beads = 200
        τ = 1.0 / (T * n_beads)

        # path = Path(n_beads, n_particles, n_dimensions, τ, m = m)
        path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)

        # Set regime
        regime = Primitive_Regime()

        """
        Set Potential Function
        """
        
        # Potential variables
        ω = 1.0
        α = 1.0
        ħ = 1.0
        
        #potential = FrohlichPotential(α,ω,ħ)
        potential = HarmonicPotential(ω)
        #potential = MexicanHatPotential(80000.0)
        #potential = ConstantPotential(10.0)

        """
        PIMC Variables
        """

        #number of steps
        n_steps = 200000

        #skipping between sampling
        equilibrium_skip = 0.5 * n_steps
        #equilibrium_skip = 0
        observables_skip = 0.001 * n_steps
        #observables_skip = 10 * n_steps

        #types of moves
        #movers = [[Bisect!],[1.0]]
        movers = [[Single!],[1.0]]
        #movers = [[Displace!,Single!],[0.2,1.0]]

        observables = [Energy]
        
        estimators = [Virial_Estimator()]
        #estimators = [Thermodynamic_Estimator()]
        #estimators = [Simple_Estimator()]
            
        """
        Run Simulation
        """

        thermalised_start!(path, potential, n_steps = 10000)
        pimc = PIMC(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true, visual=false)
        acceptance_ratio = pimc[1]
        output_observables = pimc[2]

        energy = output_observables["Energy"][string(Symbol(estimators[1]))]
        #position = output_observables["Position"][string(Symbol(estimators[1]))]

        # Comparison energy
        if typeof(potential) == HarmonicPotential
            comparison_energy = analytic_energy_harmonic(potential,1/T,ħ)
        elseif typeof(potential) == FrohlichPotential
            comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
            comparison_energy = comparison_polaron.F
        end

        println(length(energy))
        variances = jackknife(energy)

        # Post analysis
        println("acceptance ratio = ", acceptance_ratio)
        println("Mean energy = ", mean(energy))
        println("comparison_energy = ", comparison_energy)
        println("jackknife errors = ", sqrt(variances[2]))
        push!(Mean_energy_arr, mean(energy))
        push!(Comparison_energy_arr, comparison_energy)
        push!(Error_arr, sqrt(variances[2]))

        #Plots
        energyplot = plot(energy, ylabel="Energy", xlabel="x * $observables_skip steps")
        #posplot = histogram(position[:,1,1])
        #plot(posplot, energyplot, layout = (2,1), legend = false)
        #plot(posplot, xlabel="Position", ylabel="Prob Amplitude", legend = false)

        #=
        function harmonic_write(β)
            y = 1/2*ω*ħ * (1+exp(-ħ*ω*β))/(1-exp(-ħ*ω*β))
        end
        println("Analytical is ", harmonic_write(β))
        =#
 
    end

    default(fontfamily="Computer Modern",
        titlefont = (20, "Computer Modern"),
        legendfontsize = 12,
        guidefont = (18, "Computer Modern"),
        tickfont = (12, "Computer Modern"),
        linewidth=2, framestyle=:box, label=nothing, grid=false)
    #scalefontsizes()

    Beta_plot = scatter(T_arr,
                        Mean_energy_arr,
                        yaxis=:log,
                        xaxis=:log,
                        labels = "Data",
                        legend=:topleft,
                        yerror = Error_arr,
                        markerstrokewidth=0,
                        markercolor = "Blue",
                        ylabel=L"\textrm{Energy\,/\, } E",
                        xlabel=L"\textrm{Temperature\,/\,} T")
    scatter!(T_arr, Comparison_energy_arr, labels = "Theory", markerstrokewidth=0, markercolor = "Red",)
    
    title!("Energy vs Temperature")
    display(Beta_plot)
    savefig("./figs/HarmonicEnergyVSTemperature.png")
end




#=
begin

# Visualise

anim = animate_PIMC(pimc, n_particles)
gif(anim, "saved_plots/anim_output.gif", fps = 8)

end

#testing
begin
    estimators = [Thermodynamic_Estimator()]

    estimators_string = []
    for estimator in estimators
        push!(estimators_string, string(Symbol(estimator)))
    end
end
=#

