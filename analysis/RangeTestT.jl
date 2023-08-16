begin
    include("PIMC_Range.jl")

    default(fontfamily="Computer Modern",
            titlefont = (20, "Computer Modern"),
            legendfontsize = 12,
            guidefont = (18, "Computer Modern"),
            tickfont = (12, "Computer Modern"),
            linewidth=2, framestyle=:box, label=nothing, grid=false)
end

@time begin
    n_steps = 500000;
    data_set, energy_arr, error_arr = RangeTempPIMC(
                [0.1, 0.2, 0.3, 0.4, 0.5, 0.6], #Temperature
                1.0, # mass
                1.0, # ω (has to be float)
                2.0, # α (has to be float)
                1, # no of particles
                3, # number of n_dimensions
                "Primitive", # regime type
                true, # fixing beads or not
                0.03, # fixed_τ
                50, # n_beads (if τ is not fixed)
                n_steps, # No. of steps
                50000, # number of thermalisation steps
                "Single", # movers
                "Harmonic", # potential type
                "Virial" # estimators
                );
    println("-----Simulation Ended-----")
end