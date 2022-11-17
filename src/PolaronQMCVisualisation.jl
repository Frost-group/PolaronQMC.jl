# PolaronQMCVisualisation.jl
using Revise
using Plots
using LaTeXStrings


function draw_beads_3d(beads, xlims, ylims, zlims, n_particles, frame, potential, mover, T)

	p = plot(xlabel=L"x",ylabel=L"y", zlabel=L"z")
	
	for particle in 1:n_particles

		x = beads[:, particle, 1]
		y = beads[:, particle, 2]
		z = beads[:, particle, 3]
		
		x = Array(reshape(x, length(x)))
		y = Array(reshape(y, length(y)))
		z = Array(reshape(z, length(z)))

		push!(x, x[1])
		push!(y, y[1])
		push!(z, z[1])

		plot!(p, x, y, z, marker = :circle, label = "P $particle", legend = false, xlims = xlims, ylims = ylims, zlims = zlims)
		title!("PIMC Animation [Frame=$(frame)] \n $(potential) \n $(mover) \n T=$(T)", title_position=:left)

	end
	return p
end
    

function animate_PIMC(pimc, n_particles, potential, mover, T)

	xlims = [minimum([minimum(x[:, :, 1]) for x in pimc[3]]), maximum([maximum(x[:, :, 1]) for x in pimc[3]])]
	ylims = [minimum([minimum(x[:, :, 2]) for x in pimc[3]]), maximum([maximum(x[:, :, 2]) for x in pimc[3]])]
	zlims = [minimum([minimum(x[:, :, 3]) for x in pimc[3]]), maximum([maximum(x[:, :, 3]) for x in pimc[3]])]

	frame = 0
	animation = Plots.@animate for p in pimc[3]
		draw_beads_3d(p, xlims, ylims, zlims, n_particles, frame, potential, mover, T)
		frame += 1
	end

	return animation
end
