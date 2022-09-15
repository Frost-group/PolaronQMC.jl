begin
using Revise
using PolaronQMC
using BenchmarkTools
end

begin

#variables used
T = 3.0
β = 1/T
alpha_range = 1.0
thermalisation_steps = 100
steps_base = 800
fixed_τ = 0.0025
adjusted_beads = Int(floor(1/(fixed_τ*T)))

path = Path(adjusted_beads, 1, 3, fixed_τ)

α = 1.0


println("Benchmarking PIMC moves...")
println("Beads: $adjusted_beads ")

particle = rand(1:path.n_particles)

# Could set histogram limits here etc. ; see
# https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Visualizing-benchmark-results
io=IOContext(stdout)

println("Single!()")
s_bm = @benchmark Single!(path, particle, FrohlichPotential(α, 1.0, 1.0), Primitive_Regime(), path.adjusters["Single!"])
show(io,MIME("text/plain"), s_bm)

println(" ")

println("Bisect!()")
b_bm = @benchmark Bisect!(path, particle, FrohlichPotential(α, 1.0, 1.0), Primitive_Regime(), path.adjusters["Bisect!"])
show(io,MIME("text/plain"),b_bm)
end

