module PolaronQMCVisualisation

import Plots


function draw_beads_3d(path, xlims, ylims, zlims)

	p = plot()
	
	for particle in 1:path.n_particles
		x = path[:, particle, 1]
		y = path[:, particle, 2]
		z = path[:, particle, 3]
		
		x = reshape(x, length(x))
		y = reshape(y, length(y))
		z = reshape(z, length(z))

		push!(x, x[1])
		push!(y, y[1])
		push!(z, z[1])

		plot!(p, x, y, z, marker = :circle, label = "P $particle", legend = false, xlims = xlims, ylims = ylims, zlims = zlims)
	end
	return p
end

function animate_PIMC(pimc)

	xlims = [minimum([minimum(x[:, :, 1]) for x in pimc[2]]), maximum([maximum(x[:, :, 1]) for x in pimc[2]])]
	ylims = [minimum([minimum(x[:, :, 2]) for x in pimc[2]]), maximum([maximum(x[:, :, 2]) for x in pimc[2]])]
	zlims = [minimum([minimum(x[:, :, 3]) for x in pimc[2]]), maximum([maximum(x[:, :, 3]) for x in pimc[2]])]

	animation = Plots.@animate for p in pimc[2]
		draw_beads_3d(p, xlims, ylims, zlims)
	end

	return animation
end

end # sub-module
