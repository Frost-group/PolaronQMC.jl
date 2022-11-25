begin
    include("GeneralPIMC.jl")
    using Base.Threads

    default(fontfamily="Computer Modern",
            titlefont = (20, "Computer Modern"),
            legendfontsize = 12,
            guidefont = (18, "Computer Modern"),
            tickfont = (12, "Computer Modern"),
            linewidth=2, framestyle=:box, label=nothing, grid=false)
            #scalefontsizes()
end

    

begin
    #types of moves
    #movers = [[Bisect!],[1.0]]
    movers = [[Single!],[1.0]]
    #movers = [[Displace!,Single!],[0.2,1.0]]
    Potential = "Harmonic"
    #Estimator = "Thermodynamic"
    Estimator = "Virial"
    T_arr = 0.1:0.1:1.0
    #T_arr = 10 .^ range(-2, 2, length = 8)
    delta_tau = 0.2
    #T_arr = [0.2] # for a single test
    Mean_energy_arr = [0.0 for i in 1:length(T_arr)]
    Error_arr = [0.0 for i in 1:length(T_arr)]
    
    threads = true # Can turn it off by setting to false, if thread turn on cannot do plot inside function generalPIMC
    #threads = false
    if threads
        @threads for i in 1:length(T_arr)
            println("Temp is", T_arr[i])
            # println("i = $(T_arr[i]) on thread $(Threads.threadid())")
            energy, variances, accept_ratio = generalPIMC(T_arr[i], #Temperature
                                            1.0, # mass
                                            1.0, # ω (has to be float)
                                            1.0, # α (has to be float)
                                            1, # no of particles
                                            1, # number of n_dimensions
                                            Simple_Regime(), # regime type
                                            true, # fixing tau or not
                                            delta_tau, # fixed_τ
                                            200, # n_beads of tau not fixed
                                            1000000, # No. of steps
                                            100000, # number of thermalisation
                                            movers, # movers
                                            Potential, # potential type
                                            Estimator, # estimators
                                            threads # not threading
                                            )            
            #push!(Mean_energy_arr, mean(energy))
            #push!(Error_arr, sqrt(variances[2]))
            Mean_energy_arr[i] = mean(energy)
            Error_arr[i] = sqrt(variances[2])
        end
    else
        for i in 1:length(T_arr)
            println("Temp is", T_arr[i])
            # println("i = $(T_arr[i]) on thread $(Threads.threadid())")
            energy, variances, accept_ratio = generalPIMC(T_arr[i], #Temperature
                                            1.0, # mass
                                            1.0, # ω (has to be float)
                                            1.0, # α (has to be float)
                                            1, # no of particles
                                            1, # number of n_dimensions
                                            Simple_Regime(), # regime type
                                            true, # fixing tau or not
                                            delta_tau, # fixed_τ
                                            200, # n_beads of tau not fixed
                                            1000000, # No. of steps
                                            100000, # number of thermalisation
                                            movers, # movers
                                            Potential, # potential type
                                            Estimator, # estimators
                                            threads # not threading
                                            )         
            #push!(Mean_energy_arr, mean(energy))
            #push!(Error_arr, sqrt(variances[2]))
            Mean_energy_arr[i] = mean(energy)
            Error_arr[i] = sqrt(variances[2])
        end
    end
end

begin
    hbar = 1
    omega = 1
    #T_arr = 0.1:0.1:0.8
    Analytical_QHO(x) = hbar*omega/2 .+ hbar*omega*exp.(-omega*hbar ./ x)./(1 .-exp.(-hbar*omega ./ x))
    Energy_plot = plot(Analytical_QHO, 0.1, 1, labels="analytical", legend=:topleft)  
    
    scatter!(T_arr,
            Mean_energy_arr,
            labels = Estimator,
            legend=:topleft,
            xaxis=:log,
            yerror = Error_arr,
            markerstrokewidth=3,
            markercolor = "Red",
            ylabel=L"\textrm{Energy\,/\, } E",
            xlabel=L"\textrm{Temperature\,/\,} T")
    
    #scatter!(T_arr, Comparison_energy_arr, labels = "Theory", markerstrokewidth=0, markercolor = "Red",)
    title!("Energy vs Temperature")
    display(Energy_plot)
    #savefig("./figs/HarmonicEnergyVSTemperature_$(Estimator)_deltatau=$(delta_tau).png")
    #savefig("../figs/HarmonicEnergyVSTemperature_$(Estimator).png")
    #savefig("./figs/HarmonicEnergyVSTemperature_break_down.png")
end