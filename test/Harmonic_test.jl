#Testing the harmonic with different values of T

using Revise
using PolaronQMC
using Statistics
using Plots
include("../src/PolaronQMCVisualisation.jl")
using .PolaronQMCVisualisation
using PolaronMobility
using LaTeXStrings
using JLD
using Base.Threads


begin
    
    """
    Initialise System Variables
    """
    version = Int(floor(rand()*10000))
    T_arr = 0.1:0.05:0.8
    m = 1.0
    ω = 1.0
    α = 1.0
    ħ = 1.0
    n_particles = 1
    n_dimensions = 1
    start_range = 1.0


    Mean_energy_arr = [0.0 for i in 1:length(T_arr)]
    Error_arr = [0.0 for i in 1:length(T_arr)]
    Comparison_energy_arr = [0.0 for i in 1:length(T_arr)]
    data_arr = [Any[] for i in 1:length(T_arr)]

    for i in 1:length(T_arr)

        T = T_arr[i]
        β = 1 / T

        # Number of Monte-Carlo-Steps
        n_steps = 30000

        # Choose potential from "Harmonic", "Frohlich", "MexicanHat", "Constant"
        potential = "Harmonic"
        pot = potential

        # Choose Monte-Carlo Mover from "Single", "Displace", "Bisect"
        mover = "Single"

        # Choose path regime from "Simple", "Primitive"
        regime = "Primitive"
        
        # Choose energy estimator from "Simple", "Virial", "Thermodynamic"
        estimator = "Virial"
        
        # Choose Observables
        observables = ["Energy", "Position", "Correlation"]
        
        # Choose Estimators
        energy_estimators = []

        # Pick True for fixed beads or False for fixed τ
        fixed_beads = false
        if fixed_beads
            n_beads = 500
            τ = 1.0 / (T * n_beads)
        else
            # For fixed τ
            τ = 0.08
            n_beads = Int(floor(1/(τ*T)))
        end

        # Fixed observable skip or step dependant
        quick_steps = false
        if quick_steps
            equilibrium_skip = 1000
            observables_skip = 1000
        else
            #equilibrium_skip = 0.5 * n_steps #try to put as 0.5, for 0.2 for quicker testing process
            equilibrium_skip = 5000 #try to put as 0.5, for 0.2 for quicker testing process
            observables_skip = 0.005 * n_steps
        end

        # Initate path
        path = Path(n_beads, n_particles, n_dimensions, τ, m=m)

        # Set regime
        if regime == "Primitive"
            regime = PrimitiveRegime()
        elseif regime == "Simple"
            regime = SimpleRegime()
        elseif regime == "LBRegime"
            regime = LBRegime()
        elseif regime == "BoundRegime"
            regime = BoundRegime()
        else
            println("Invalid Regime: ", regime)
        end
    
        # Set Potential
        if potential == "Frohlich"
            potential = FrohlichPotential(α,ω,ħ)
        elseif potential == "Harmonic"
            potential = HarmonicPotential(ω)
        elseif potential == "MexicanHat"
            potential = MexicanHatPotential(80.0)
        elseif potential == "Contsant"
            potential = ConstantPotential(10.0)
        else
            println("Invalid Potential: ")
        end

        # Set Monte-Carlo Mover
        if mover == "Single"
            mover = SingleMover(path)
        elseif mover == "Displace"
            mover = DisplaceMover(path)
        elseif mover == "Bisect"
            mover = BisectMover(path)
        else
            println("Invalid mover")
        end

        # Set Estimator
        if estimator == "Virial"
            estimators = [VirialEstimator()]
        elseif estimator == "Thermodynamic"
            estimators = [ThermodynamicEstimator()]
        elseif estimator == "Simple"
            estimators = [SimpleEstimator()]
        else
            println("Invalid Estimator: ", estimator)
        end

        """
        Run Simulation
        """

        # thermalised_start!(path, potential, n_steps = 100000)
        data = PIMC(n_steps, equilibrium_skip, observables_skip, path, mover, estimators, potential, regime, observables, adjust=true)
        
        # Store outputs
        energies = data["Energy:$(estimator)"]
        positions = data["Position:p1d1"]
        correlations = data["Correlation:$(estimator)"]
            
        if typeof(potential) == HarmonicPotential
            comparison_energy = analyticEnergyHarmonic(potential,β,ħ,n_dimensions) * n_particles
        elseif typeof(potential) == FrohlichPotential
            comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true) * n_particles
            comparison_energy = comparison_polaron.F
        end

        # Post analysis
        variances = jackknife(energies)
        jacknife_errors = sqrt(variances[2])
        mean_energy = mean(energies)
        corr_mean = mean(correlations)
        corr_std = std(correlations)

        Mean_energy_arr[i] = mean_energy
        Error_arr[i] = jacknife_errors
        Comparison_energy_arr[i] = comparison_energy
        data_arr[i] = data
    end

    save("data_arr/Harmonic/$(string(typeof(potential)))_T$(T)_nsteps$(n_steps)_v$(version)_beads$(n_beads).jld", "data_arr", data_arr, "energies", energies, "comparison_energy", comparison_energy, "correlations", correlations, "jacknife_errors", jacknife_errors, 
            "equilibrium_skip", equilibrium_skip, "observables_skip", observables_skip)

    default(fontfamily="Computer Modern",
        titlefont = (20, "Computer Modern"),
        legendfontsize = 12,
        guidefont = (18, "Computer Modern"),
        tickfont = (12, "Computer Modern"),
        linewidth=2, framestyle=:box, label=nothing, grid=false)
    #scalefontsizes()

    Beta_plot = scatter(T_arr,
                        Mean_energy_arr,
                        labels = "data",
                        legend=:topleft,
                        yerror = Error_arr,
                        markerstrokewidth=3,
                        markercolor = "Blue",
                        ylabel=L"\textrm{Energy\,/\, } E",
                        xlabel=L"\textrm{Temperature\,/\,} T")
    scatter!(T_arr, Comparison_energy_arr, labels = "Theory", markerstrokewidth=0, markercolor = "Red",)
    
    title!("Energy vs Temperature")
    display(Beta_plot)
    #savefig("./figs/HarmonicEnergyVSTemperature.png")
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

