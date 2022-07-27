T = 1.0
λ = 0.5
n_beads = 100
τ = 1.0 / (T * n_beads)
n_particles = 1
path = PolaronQMC.Path(n_beads, n_particles, τ = τ, λ = λ)

n_steps = 10_000
pimc = PolaronQMC.PIMC(n_steps, path, [Single!, Displace!], [Energy, Correlation], HarmonicPotential(1.0))

#sampled_energy = mean(pimc[3]["Energy"])

