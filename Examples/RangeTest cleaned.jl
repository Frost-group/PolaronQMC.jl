begin
    using Revise
    using PolaronQMC
    using Statistics
    using Plots
    using PolaronMobility
end

begin
    T = 3.0
    β = 1/T
    alpha_range = 1.0:0.5:3.0
    thermalisation_steps = 5000
    steps_base = 800
    fixed_τ = 0.0002
    quickrun_results = quickrun_frohlich(T, alpha_range, fixed_τ, thermalisation_steps, steps_base, thermalised = true)

end
begin
#Plotting and comparison


comparison_polaron = make_polaron(alpha_range, [T], [0.0]; ω = 1.0, rtol = 1e-4, verbose = false, threads = true)
comparison_range_L = comparison_polaron.F

addon = 1.5 * coth(β) #phonon energy addon
plot(alpha_range, comparison_range_L, label="Comparison Energy", xlabel="α", ylabel="Energy", linestyle=:dash,linecolor=:red, linewidth = 1.5, legend=:bottomleft)
scatter!(alpha_range, quickrun_results[1] .+ addon , yerr = quickrun_results[2], label="PIMC Results")
end


