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
    version = Int(floor(rand() * 10000))
    # Path variables
    n_steps = 150
    m = 1.0
    α_arr = 1.0:1.0:6.0
    Mean_energy_arr = [0.0 for i in 1:length(α_arr)]
    Error_arr = [0.0 for i in 1:length(α_arr)]
    n_particles = 1
    n_dimensions = 3

    @threads for i in 1:length(α_arr)
        
        T = 0.1
        start_range = 1.0
        β = 1 / T

        # For fixed τ
        fixed_τ = 0.05
        adjusted_beads = Int(floor(1/(fixed_τ*T)))

        # For fixed number of beads
        n_beads = 1500
        τ = 1.0 / (T * n_beads)

        path = Path(n_beads, n_particles, n_dimensions, τ)
        #path = Path(adjusted_beads, n_particles, n_dimensions, fixed_τ)

        # Set regime
        regime = Primitive_Regime()

        """
        Set Potential Function
        """
        
        # Potential variables
        ω = 1.0
        α = α_arr[i]
        ħ = 1.0
        
        potential = FrohlichPotential(α,ω,ħ)
        #potential = HarmonicPotential(ω)
        #potential = MexicanHatPotential(80.0)
        #potential = ConstantPotential(10.0)

        """
        PIMC Variables
        """

        # number of steps



        #skipping between sampling
        #equilibrium_skip = 0.1 * n_steps
        equilibrium_skip = 10
        #observables_skip = 0.001 * n_steps
        observables_skip = 5

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
    α_arr = 1.0:1.0:6.0
    Estimator = "Virial"
    #T_arr = 0.1:0.1:0.8
    #Analytical_QHO(x) = ħ*ω/2 .+ ħ*ω*exp.(-ω*ħ ./ x)./(1 .-exp.(-ħ*ω ./ x))
    
    #Analytical_QHO(x) = n_dimensions/2*ħ*ω * (1 .+ exp.(-ħ*ω ./ x))./(1 .- exp.(-ħ*ω ./ x))
    VMC = [-0.4428, -1.501, -2.623, -3.80, -5.06, -6.42]
    Energy_plot = plot(α_arr, VMC, labels="Analytical-Froh", legend=:topleft)  
    
    scatter!(α_arr,
            Mean_energy_arr,
            labels = Estimator,
            legend=:bottomleft,
            #xaxis=:log,
            yerror = Error_arr,
            markerstrokewidth=3,
            markercolor = "Red",
            ylabel=L"\textrm{Energy\,/\, } E",
            xlabel=L"\textrm{α}")
    
    #scatter!(T_arr, Comparison_energy_arr, labels = "Theory", markerstrokewidth=0, markercolor = "Red",)
    title!("Energy vs α")
    display(Energy_plot)
    savefig("./figs/Froh/3D-Frohlich_EVST_$(Estimator)_deltatau=0.05_nsteps=$(n_steps)_v$(version).png")
end

begin
    # Plots
    energy_plot = plot(energies, ylabel="Energy", xlab = "Sweeps / $observables_skip\$ n\$")
    hline!([comparison_energy], linestyle=:dash)
    energy_hist = histogram(energies, ylab="Frequencies", xlab="Energy")
    acceptance_rate_plot = plot(acceptance_rates[Int(length(acceptance_rates)*0.9):end], xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Acceptance\, Rate\, /\, } r", dpi=600)
    shift_width_plot = plot(shift_widths, xlab = L"\mathrm{Sweeps\, /\, } n", ylab=L"\mathrm{Shift\, Width\, /\, } \Delta x", dpi=600)
    acceptance_shift_plot = scatter(acceptance_rates, shift_widths, xlab=L"\mathrm{Acceptance\, Rate\, /\, } r", ylab=L"\mathrm{Shift\, Width\, /\, } \Delta x", dpi=600)
    acceptance_rate_hist = histogram(acceptance_rates, ylab="Frequency", xlab=L"\mathrm{Acceptance\, Rate\, /\, } r")
    #shift_width_hist = histogram(shift_widths, ylab="Frequency", xlab=L"\mathrm{Shift\, Width\, /\, } \Delta x")
    #posplot = histogram(position[:,1,1])
    #plot(posplot, energyplot, layout = (2,1), legend = false)
    #plot(posplot, xlabel="Position", ylabel="Prob Amplitude", legend = false)

    display(energy_hist)
    display(energy_plot)
    display(acceptance_rate_plot)
    display(shift_width_plot)
    display(acceptance_shift_plot)
    #display(shift_width_hist)
    display(acceptance_rate_hist)

    #savefig(energy_hist, "saved_plots/energy_hist.png") 
    #savefig(energy_plot, "saved_plots/energy_plot.png")    
    #savefig(acceptance_rate_plot, "saved_plots/acceptance_rate_convergence_02_12.png")
    #savefig(shift_width_plot, "saved_plots/shift_width_convergence_02_12.png")
    #savefig(acceptance_shift_plot, "saved_plots/acceptance_shift_02_12.png")
    #savefig(shift_width_hist, "saved_plots/shift_width_hist_02_12.png")
    #savefig(acceptance_rate_hist, "saved_plots/acceptance_rate_shift_02_12.png")

    
    posplot = histogram(position1[:,1,1])
    display(posplot)
end