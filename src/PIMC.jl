# PIMC.jl

function relabel_beads!(path::Path)
	rand_slice = rand(1:path.n_beads)
	slices = vcat(rand_slice:path.n_beads, 1:rand_slice-1)
	path.beads = path.beads[slices, :, :]
end

function PIMC(n_steps::Int, path::Path, movers, observables, potentials::Union{Potential, Array{Potential}})
	
	observable_skip = 0.01 * n_steps
	equilibrium_skip = 0.1 * n_steps
	
	n_accepted = Dict(string(Symbol(mover!)) => 0 for mover! in movers)
	# observable_traces = Dict(string(Symbol(observable)) => [] for observable in observables)
	
	for observable in observables
		@eval $(Symbol(observable, "_statistics")) = Series(Trace(Mean()), Trace(AutoCov(0)))
	end

	ymax = maximum(path.beads[:, :, 1])
	ymin = minimum(path.beads[:, :, 1])

	path_trace = []
	for step in 1:n_steps

		for mover! in movers
			for particle in rand(1:path.n_particles, path.n_particles)
				n_accepted[string(Symbol(mover!))] += mover!(path, particle, potentials)
			end
		end

		if mod(step, observable_skip) == 0 && step > equilibrium_skip
			observable_plots = []
			for observable in observables
				@eval trace_symbol = $(Symbol(observable, "_statistics"))
				trace_value = observable(path, potentials)
				# append!(observable_traces[string(Symbol(observable))], [trace_value])
				fit!(trace_symbol, trace_value)
				append!(observable_plots, [plot(trace_symbol, minorgrid = true, linestyle = :solid, title = "$observable", legendfontsize=5, markersize = 2.5)])
			end
			append!(path_trace, [path.beads])
			ymax = ymax < maximum(path.beads[:, :, 1]) ? maximum(path.beads[:, :, 1]) : ymax
			ymin = ymin > minimum(path.beads[:, :, 1]) ? minimum(path.beads[:, :, 1]) : ymin
			path_plot = plot(1:path.n_beads, path.beads[:, :, 1], grid = false, ylabel = "x position", xlabel = "imaginary time", ylims = [ymin, ymax], title = "Path", marker = :circle, legend = false)
			plot(observable_plots..., path_plot, size = (500, 500),  layout = @layout [a;b;c])
			gui()
		end
	end

	acceptance_ratio = Dict(string(Symbol(mover!)) => 1.0 * n_accepted[string(Symbol(mover!))] / (n_steps * path.n_particles) for mover! in movers)
	
	returns = [acceptance_ratio, path_trace]

	for observable in observables
		@eval trace_symbol = $(Symbol(observable, "_statistics")) 
		append!(returns, [trace_symbol])
	end

	return tuple(returns...)
end

