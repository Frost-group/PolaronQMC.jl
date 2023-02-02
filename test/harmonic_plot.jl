using Revise
using PolaronQMC
using Statistics
using Plots
using PolaronMobility
using LaTeXStrings
using Base.Threads

default(fontfamily="Computer Modern",
            titlefont = (16, "Computer Modern"),
            guidefont = (18, "Computer Modern"),
            tickfont = (12, "Computer Modern"),
            legendfontsize = 12,
            linewidth=2, framestyle=:box, label=nothing, grid=true)

@time begin
    

    """
    Initialise System Variables
    """
    version = 3
    # Path variables
    n_steps = 20000
    m = 1.0
    T_arr = 0.1:0.1:0.9
    Mean_energy_arr = [0.0 for i in 1:length(T_arr)]
    Error_arr = [0.0 for i in 1:length(T_arr)]
    n_particles = 1
    n_dimensions = 3

    @threads for i in 1:length(T_arr)
        
        T = T_arr[i]
        start_range = 1.0
        β = 1 / T

        # For fixed τ
        fixed_τ = 0.05
        adjusted_beads = Int(floor(1/(fixed_τ*T)))

        # For fixed number of beads
        n_beads = 50
        τ = 1.0 / (T * n_beads)

        #path = Path(n_beads, n_particles, n_dimensions, τ)
        path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)

        # Set regime
        regime = Primitive_Regime()

        """
        Set Potential Function
        """
        
        # Potential variables
        ω = 1.0
        α = 2.0
        ħ = 1.0
        
        #potential = FrohlichPotential(α,ω,ħ)
        potential = HarmonicPotential(ω)
        #potential = MexicanHatPotential(80.0)
        #potential = ConstantPotential(10.0)

        """
        PIMC Variables
        """

        # number of steps



        #skipping between sampling
        #equilibrium_skip = 0.1 * n_steps
        equilibrium_skip = 100
        #observables_skip = 0.001 * n_steps
        observables_skip = 50

        # types of moves
        movers = Dict("Bisect!" => [1.0])
        #movers = Dict("Single!" => [1.0])
        #movers = Dict("Displace!" => [1.0])
        #movers = Dict("Single!" => [1.0], "Bisect!" => [0.2])

        observables = [Energy, Position]
        
        estimators = [Virial_Estimator()]
        #estimators = [Thermodynamic_Estimator()]
        #estimators = [Simple_Estimator()]
    
        #initial_pos = Array(path.beads)

        """
        Run Simulation
        """

        # thermalised_start!(path, potential, n_steps = 100000)
        pimc = PIMC(n_steps, equilibrium_skip, observables_skip, path, movers, observables, estimators, potential, regime, adjust=true, visual=true)
        
        # Store outputs
        adjuster_stats = pimc[1]
        output_observables = pimc[2]

        mover = "Bisect!"
        #mover = "Single!"

        energies = output_observables["Energy"][string(Symbol(estimators[1]))]
        acceptance_rates = adjuster_stats[mover]["Acceptance Rate"]
        shift_widths = adjuster_stats[mover]["Shift Width"]
        position1 = output_observables["Position"][string(Symbol(estimators[1]))]

        # Comparison energy
        if typeof(potential) == HarmonicPotential
            comparison_energy = analytic_energy_harmonic(potential,β,ħ,n_dimensions)
        elseif typeof(potential) == FrohlichPotential
            comparison_polaron = make_polaron([α], [T], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true)
            comparison_energy = comparison_polaron.F
        end

        
        # Post analysis
        variances = jackknife(energies)
        jacknife_errors = sqrt(variances[2])
        mean_energy = mean(energies)
        last_acceptance_rate = last(acceptance_rates)
        mean_acceptance_rate = mean(acceptance_rates)
        std_acceptance_rate = std(acceptance_rates)

        Mean_energy_arr[i] = mean_energy
        Error_arr[i] = jacknife_errors

        # Output measurements and statistics
        #println("Number of Beads: ", adjusted_beads)
        println("Number of Beads: ", n_beads)
        println("Mean Energy: ", mean_energy)
        println("Comparison Energy: ", comparison_energy)
        println("jackknife errors: ", jacknife_errors)
        println("Final Acceptance Rate: ", last_acceptance_rate)
        println("Mean Acceptance Rate: ", mean_acceptance_rate, " +/- ", std_acceptance_rate)
    
    
    end
end

begin
    ħ = 1.0
    ω = 1.0
    Estimator = "Virial"
    #T_arr = 0.1:0.1:0.8
    #Analytical_QHO(x) = ħ*ω/2 .+ ħ*ω*exp.(-ω*ħ ./ x)./(1 .-exp.(-ħ*ω ./ x))
    
    Analytical_QHO(x) = n_dimensions/2*ħ*ω * (1 .+ exp.(-ħ*ω ./ x))./(1 .- exp.(-ħ*ω ./ x))
    
    Energy_plot = plot(Analytical_QHO, 0.1, 1, labels="analytical", legend=:topleft)  
    
    scatter!(T_arr,
            Mean_energy_arr,
            labels = Estimator,
            legend=:topleft,
            #xaxis=:log,
            yerror = Error_arr,
            markerstrokewidth=3,
            markercolor = "Red",
            ylabel=L"\textrm{Energy\,/\, } E",
            xlabel=L"\textrm{Temperature\,/\,} T")
    
    #scatter!(T_arr, Comparison_energy_arr, labels = "Theory", markerstrokewidth=0, markercolor = "Red",)
    title!("Energy vs Temperature")
    display(Energy_plot)
    savefig("./figs/QHO/3D-QHO_EVST_$(Estimator)_deltatau=0.05_nsteps=$(n_steps)_v$(version).png")
end

