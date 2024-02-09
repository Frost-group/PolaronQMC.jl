using PolaronQMC
using Plots
using PolaronMobility
using LaTeXStrings

# Comparison data for Frohlich

function WeakPTFrohlich(alpha)
    return -1 * alpha - 0.0159196 * alpha^2
end

function strongPTFrohlich(alpha)
    return -2.84 - 0.1085 * alpha^2 #Pertubation theory
end

α_range = 1.0:0.1:12.0

DiagMC_α = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
DiagMC_α_E =
    [-1.08, -2.13, -3.2, -4.36, -5.59, -6.91, -8.42, -9.98, -11.77, -13.8, -16.05, -18.52]

Var_data = []
variational_run = feynmanvw.(2.0, 1.0, collect(α_range), 1.0, Inf)
for i = 1:length(variational_run)
    push!(Var_data, variational_run[i][3])
end

Frohlich_plot = plot(
    α_range,
    WeakPTFrohlich.(α_range),
    labels = "WeakPT",
    linestyle = :dash,
    ylabel = L"Energy \ E(α) \ (T=0.1)",
    xlabel = L"α",
    dpi = 300,
)
plot!(
    α_range,
    strongPTFrohlich.(α_range),
    labels = "StrongPT",
    linestyle = :dash,
    dpi = 300,
)
plot!(α_range, Var_data, label = "Var")
scatter!(DiagMC_α, DiagMC_α_E, label = "DiagMC", markersize = :3)
display(Frohlich_plot)
