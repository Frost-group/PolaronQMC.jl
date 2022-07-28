# PIMC.jl

function relabel_beads!(path::Path)
	rand_slice = rand(1:path.n_beads)
	slices = vcat(rand_slice:path.n_beads, 1:rand_slice-1)
	path.beads = path.beads[slices, :, :]
end

function PIMC(n_steps::Int, path::Path, movers, observables, potentials::Union{Potential, Array{Potential}}; n_particles = 1, observable_skip = 10, equilibrium_skip = 1000, visual = true)
		
	n_accepted = Dict(string(Symbol(mover!)) => 0 for mover! in movers)
	
	for observable in observables
		@eval $(Symbol(observable, "_statistics")) = Series(Trace(Mean()), Trace(AutoCov(0)))
	end

	if visual
		ymax = maximum(path.beads[:, :, 1])
		ymin = minimum(path.beads[:, :, 1])
	end

	path_trace = []
	for step in 1:n_steps

		for mover! in movers
			for particle in rand(1:path.n_particles, n_particles)
				n_accepted[string(Symbol(mover!))] += mover!(path, particle, potentials)
			end
		end

		if mod(step, observable_skip) == 0 && step > equilibrium_skip
			observable_plots = []

			for observable in observables
				@eval trace_symbol = $(Symbol(observable, "_statistics"))
				trace_value = observable(path, potentials)
				fit!(trace_symbol, trace_value)

				if visual
					append!(observable_plots, [plot(trace_symbol, minorgrid = true, linestyle = :solid, title = "$observable", legendfontsize=5, markersize = 2.5)])
				end
			end

			append!(path_trace, [path.beads])
			if visual
				ymax = ymax < maximum(path.beads[:, :, 1]) ? maximum(path.beads[:, :, 1]) : ymax
				ymin = ymin > minimum(path.beads[:, :, 1]) ? minimum(path.beads[:, :, 1]) : ymin
				path_plot = plot(1:path.n_beads, path.beads[:, :, 1], grid = false, ylabel = "x position", xlabel = "imaginary time", ylims = [ymin, ymax], title = "Path", marker = :circle, legend = false)
				plot(observable_plots..., path_plot, size = (1000, 500))
				gui()
			end
		end
	end

	acceptance_ratio = Dict(string(Symbol(mover!)) => 1.0 * n_accepted[string(Symbol(mover!))] / (n_steps * n_particles) for mover! in movers)
	
	returns = [acceptance_ratio, path_trace]

	for observable in observables
		@eval trace_symbol = $(Symbol(observable, "_statistics")) 
		append!(returns, [trace_symbol])
	end

	return tuple(returns...)
end

