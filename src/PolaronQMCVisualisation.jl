module PolaronQMCVisualisation

using Revise
using Plots

export animate_PIMC, draw_beads_3d


function draw_beads_3d(beads, xlims, ylims, zlims, n_particles)

	p = plot()
	
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
	end
	return p
end
    
function animate_PIMC(pimc, n_particles)

	xlims = [minimum([minimum(x[:, :, 1]) for x in pimc[3]]), maximum([maximum(x[:, :, 1]) for x in pimc[3]])]
	ylims = [minimum([minimum(x[:, :, 2]) for x in pimc[3]]), maximum([maximum(x[:, :, 2]) for x in pimc[3]])]
	zlims = [minimum([minimum(x[:, :, 3]) for x in pimc[3]]), maximum([maximum(x[:, :, 3]) for x in pimc[3]])]

	animation = Plots.@animate for p in pimc[3]
		draw_beads_3d(p, xlims, ylims, zlims, n_particles)
	end

	return animation
end


end # sub-module
