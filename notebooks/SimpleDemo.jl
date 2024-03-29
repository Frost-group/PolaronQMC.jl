### A Pluto.jl notebook ###
# v0.19.5

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

# ╔═╡ 981f58c8-fa4b-476c-a7cb-907b6fde5c4f
using BenchmarkTools

# ╔═╡ 8cae5b2b-ce15-45bd-a1ac-86fd116885b5
using Plots, Pluto

# ╔═╡ 8a6d45d3-9a97-4f0e-b484-de073fc92091
begin
    T = 0.5
    λ = 0.5
    n_beads = 8
    τ = 1.0 / (T * n_beads)
    n_particles = 1
    path = PolaronQMC.Path(n_beads, n_particles, n_dimensions = 2, τ = τ, λ = λ)

    n_steps = 1_000
    p_beads = copy(path.beads)
    pimc = PolaronQMC.PIMC(
        n_steps,
        path,
        [Single!, Displace!],
        [Energy, Correlation],
        MexicanHatPotential(5.0),
    )
end

# ╔═╡ ac76e1d5-4acc-453d-a02d-e1e3cf4feff2
sampled_energy = mean(pimc[3]["Energy"])

# ╔═╡ 57874a17-d406-46cd-ba6a-3aae950a7e01
correlation = [c[i] for i = 1:n_beads, c in (pimc[3]["Correlation"])]

# ╔═╡ ffdcb683-ab78-4a08-9ad2-1adeb6991afb
sampled_correlation = [mean(correlation[:, i]) for i = 1:n_beads]

# ╔═╡ 7ad85b78-0396-48f4-adc9-00fb5344d114
md"""
# Benchmarking internal functions reactively
"""

# ╔═╡ 16c2cfe7-c04a-45fe-ba35-ad917b3ef888
@benchmark PolaronQMC.PIMC(
    100,
    path,
    [Single!, Displace!],
    [Energy, Correlation],
    HarmonicPotential(1.0),
)

# ╔═╡ 96957c39-cb07-4894-a1d3-377477fe3ae4
md"""
# Animation
"""

# ╔═╡ 3852558b-3437-40f3-9fc9-75555b467ba6
function draw_beads_2d(path, xlims, ylims)

    p = plot()

    for particle = 1:size(path)[2]
        x = path[:, particle, 1]
        y = path[:, particle, 2]

        x = reshape(x, length(x))
        y = reshape(y, length(y))

        push!(x, x[1])
        push!(y, y[1])

        plot!(
            p,
            x,
            y,
            marker = :circle,
            label = "P $particle",
            legend = false,
            xlims = xlims,
            ylims = ylims,
        )
    end
    return p
end

# ╔═╡ 5295edba-c08d-4a26-be56-bc24ae0345aa
function draw_beads_3d(path, xlims, ylims, zlims)

    p = plot()

    for particle = 1:size(path)[2]
        x = path[:, particle, 1]
        y = path[:, particle, 2]
        z = path[:, particle, 3]

        x = reshape(x, length(x))
        y = reshape(y, length(y))
        z = reshape(z, length(z))

        push!(x, x[1])
        push!(y, y[1])
        push!(z, z[1])

        plot!(
            p,
            x,
            y,
            z,
            marker = :circle,
            label = "P $particle",
            legend = false,
            xlims = xlims,
            ylims = ylims,
            zlims = zlims,
        )
    end
    return p
end

# ╔═╡ 548bdf3b-5089-471d-8e18-192e3e2172fc
function animate_PIMC2d(pimc)

    xlims = [
        minimum([minimum(x[:, :, 1]) for x in pimc[2]]),
        maximum([maximum(x[:, :, 1]) for x in pimc[2]]),
    ]
    ylims = [
        minimum([minimum(x[:, :, 2]) for x in pimc[2]]),
        maximum([maximum(x[:, :, 2]) for x in pimc[2]]),
    ]

    animation = @animate for p in pimc[2]
        draw_beads_2d(p, xlims, ylims)
    end

    return animation
end


# ╔═╡ 05153c0f-6359-46d7-a928-f9aa50ea8113
function animate_PIMC3d(pimc)

    xlims = [
        minimum([minimum(x[:, :, 1]) for x in pimc[2]]),
        maximum([maximum(x[:, :, 1]) for x in pimc[2]]),
    ]
    ylims = [
        minimum([minimum(x[:, :, 2]) for x in pimc[2]]),
        maximum([maximum(x[:, :, 2]) for x in pimc[2]]),
    ]
    zlims = [
        minimum([minimum(x[:, :, 3]) for x in pimc[2]]),
        maximum([maximum(x[:, :, 3]) for x in pimc[2]]),
    ]

    animation = @animate for p in pimc[2]
        draw_beads_3d(p, xlims, ylims, zlims)
    end

    return animation
end


# ╔═╡ 0b18d25b-0a59-41a2-8141-00e565ae9d6d


# ╔═╡ f69336a9-9d4a-4bb6-a665-cfa1ccb50ad0
anim = animate_PIMC2d(pimc)

# ╔═╡ ab39844d-c9c0-4ef3-ac20-c091d8654626
gif(anim, "anim_fps15.gif", fps = 15)

# ╔═╡ e229663c-964a-4cbf-afad-fd34b8e18b59
pimc[2][1][:, :, 1]

# ╔═╡ 6f2921af-ddd6-43d7-ae6b-fb0a75a8c172
@code_warntype Single!(path, 1, HarmonicPotential(1.0))

# ╔═╡ Cell order:
# ╠═3f6960d4-d770-11ec-35d7-fd7da2e01cc8
# ╠═8a6d45d3-9a97-4f0e-b484-de073fc92091
# ╠═ac76e1d5-4acc-453d-a02d-e1e3cf4feff2
# ╠═57874a17-d406-46cd-ba6a-3aae950a7e01
# ╠═ffdcb683-ab78-4a08-9ad2-1adeb6991afb
# ╟─7ad85b78-0396-48f4-adc9-00fb5344d114
# ╠═981f58c8-fa4b-476c-a7cb-907b6fde5c4f
# ╠═16c2cfe7-c04a-45fe-ba35-ad917b3ef888
# ╟─96957c39-cb07-4894-a1d3-377477fe3ae4
# ╠═8cae5b2b-ce15-45bd-a1ac-86fd116885b5
# ╠═3852558b-3437-40f3-9fc9-75555b467ba6
# ╠═5295edba-c08d-4a26-be56-bc24ae0345aa
# ╠═548bdf3b-5089-471d-8e18-192e3e2172fc
# ╠═05153c0f-6359-46d7-a928-f9aa50ea8113
# ╠═0b18d25b-0a59-41a2-8141-00e565ae9d6d
# ╠═f69336a9-9d4a-4bb6-a665-cfa1ccb50ad0
# ╠═ab39844d-c9c0-4ef3-ac20-c091d8654626
# ╠═e229663c-964a-4cbf-afad-fd34b8e18b59
# ╠═6f2921af-ddd6-43d7-ae6b-fb0a75a8c172
