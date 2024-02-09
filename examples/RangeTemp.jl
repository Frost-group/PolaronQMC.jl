using PolaronQMC
using Plots
begin
    default(
        fontfamily = "Computer Modern",
        titlefont = (20, "Computer Modern"),
        legendfontsize = 12,
        guidefont = (18, "Computer Modern"),
        tickfont = (12, "Computer Modern"),
        linewidth = 2,
        framestyle = :box,
        label = nothing,
        grid = false,
    )
end

@time begin
    n_steps = 100000
    energy_arr, error_arr, data_set = RangeTempPIMC(
        [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        1.0,
        2.0,
        3,
        n_steps,
        pot = "Harmonic",
    )
    println("-----Simulation Ended-----")
end
