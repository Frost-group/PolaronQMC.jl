### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 3f6960d4-d770-11ec-35d7-fd7da2e01cc8
begin
	using Pkg
	Pkg.activate("../")
	using PolaronQMC
	using Statistics
	#PolaronQMCVisualisation
end

# ╔═╡ 8a6d45d3-9a97-4f0e-b484-de073fc92091
begin
	T = 1.0
	λ = 0.5
	n_beads = 100
	τ = 1.0 / (T * n_beads)
	n_particles = 1
	path = PolaronQMC.Path(n_beads, n_particles, τ = τ, λ = λ)

	n_steps = 10_000
	p_beads = copy(path.beads)
	pimc = PolaronQMC.PIMC(n_steps, path, [Single, Displace], [Energy, Correlation], HarmonicPotential(1.0))
end

# ╔═╡ ac76e1d5-4acc-453d-a02d-e1e3cf4feff2
sampled_energy = mean(pimc[3]["Energy"])

# ╔═╡ 57874a17-d406-46cd-ba6a-3aae950a7e01
correlation = [c[i] for i in 1:n_beads, c in (pimc[3]["Correlation"])]

# ╔═╡ ffdcb683-ab78-4a08-9ad2-1adeb6991afb
sampled_correlation = [mean(correlation[:, i]) for i in 1:n_beads]

# ╔═╡ Cell order:
# ╠═3f6960d4-d770-11ec-35d7-fd7da2e01cc8
# ╠═8a6d45d3-9a97-4f0e-b484-de073fc92091
# ╠═ac76e1d5-4acc-453d-a02d-e1e3cf4feff2
# ╠═57874a17-d406-46cd-ba6a-3aae950a7e01
# ╠═ffdcb683-ab78-4a08-9ad2-1adeb6991afb
