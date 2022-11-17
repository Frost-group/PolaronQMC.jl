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
    T_arr = 0.1:0.1:0.8
    Mean_energy_arr = [0.0 for i in 1:length(T_arr)]
    Error_arr = [0.0 for i in 1:length(T_arr)]
    threads = true # Can turn it off by setting to false, if thread turn on cannot do plot inside function generalPIMC
    #Mean_energy_arr = []
    #Error_arr = []
    if threads
        @threads for i in 1:length(T_arr)
            println("Temp is", T_arr[i])
            # println("i = $(T_arr[i]) on thread $(Threads.threadid())")
            energy, variances = generalPIMC(T_arr[i], #Temperature
                        1, # mass
                        1, # no of particles
                        1, # number of n_dimensions
                        Simple_Regime(), # regime type
                        true, # fixing tau or not
                        0.2, # fixed_τ
                        200, #n_beads if tau not fixed
                        10000000, # No. of steps
                        100000, # number of thermalisation
                        movers, # moving types
                        threads # whether thread or not
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
            energy, variances = generalPIMC(T_arr[i], #Temperature
                        1, # mass
                        1, # no of particles
                        1, # number of n_dimensions
                        Simple_Regime(), # regime type
                        true, # fixing tau or not
                        0.2, # fixed_τ
                        200, #n_beads if tau not fixed
                        10000000, # No. of steps
                        100000, # number of thermalisation
                        movers
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
    T_arr = 0.1:0.1:0.8
    Analytical_QHO(x) = hbar*omega/2 .+ hbar*omega*exp.(-omega*hbar ./ x)./(1 .-exp.(-hbar*omega ./ x))
    Energy_plot = plot(Analytical_QHO, 0.1, 1, labels="analytical", legend=:topleft)  
    
    scatter!(T_arr,
            Mean_energy_arr,
            labels = "data",
            legend=:topleft,
            yerror = Error_arr,
            markerstrokewidth=3,
            markercolor = "Red",
            ylabel=L"\textrm{Energy\,/\, } E",
            xlabel=L"\textrm{Temperature\,/\,} T")
    
    #scatter!(T_arr, Comparison_energy_arr, labels = "Theory", markerstrokewidth=0, markercolor = "Red",)
    title!("Energy vs Temperature")
    display(Energy_plot)
    #savefig("../figs/HarmonicEnergyVSTemperature.png")
end