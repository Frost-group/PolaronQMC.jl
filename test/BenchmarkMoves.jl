using BenchmarkTools

T = 1.0
λ = 0.5
n_beads = 100
τ = 1.0 / (T * n_beads)
n_particles = 1
path = PolaronQMC.Path(n_beads, n_particles, τ = τ, λ = λ)

println("Benchmarking PIMC moves...")
println("Particles: $n_particles Beads: $n_beads ")

# Could set histogram limits here etc. ; see
# https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Visualizing-benchmark-results
io=IOContext(stdout)

println("Single!()")
b= @benchmark Single!(path,1,HarmonicPotential(1.0)) 
show(io,MIME("text/plain"),b)

println("Displace!()")
b= @benchmark Displace!(path,1,HarmonicPotential(1.0))
show(io,MIME("text/plain"),b)

#sampled_energy = mean(pimc[3]["Energy"])

