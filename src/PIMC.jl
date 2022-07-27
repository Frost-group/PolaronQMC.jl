# PIMC.jl



function PIMC(n_steps::Int, equilibrium_skip, observable_skip, path::Path, movers, observables, estimators, potential::Potential, regime::Regime; adjust = true, visual = false)
	#setting up storage of output

	# information about the running of the system
		system_stats = Dict()
		system_stats["acceptance_array"] = Dict()
		system_stats["attempted_array"] = Dict()

	#acceptance and attempt arrays for movers

		for mover in movers[1]
			system_stats["acceptance_array"][string(Symbol(mover))] = 0
			system_stats["attempted_array"][string(Symbol(mover))] = 0
		end

	#output arrays for different estimators of observables
		output_observables = Dict()
		#observables_counter = 0 #used in weighting of the rolling averages of observables
		#generating lists for output
		for observable in observables
			output_observables[string(observable)] = Dict() 
			for estimator in estimators
				output_observables[string(observable)][string(Symbol(estimator))] = []
			end
		end

		

	
	for step in 1:n_steps
		if step == 0.5*n_steps
			println("50% complete")
		end

		#updating n_accepted, moving beads, and changing shift width if necessary
		for particle in rand(1:path.n_particles, path.n_particles)
			for mover_index in 1:length(movers[1])
				adjuster = path.adjusters[string(Symbol(movers[1][mover_index]))]
				if movers[2][mover_index] > rand()
					system_stats["attempted_array"][string(Symbol(movers[1][mover_index]))] += 1
					system_stats["acceptance_array"][string(Symbol(movers[1][mover_index]))] += movers[1][mover_index](path, particle, potential, regime, adjuster)
					

				end
			end

			if adjust #changing shift width automatically
				for adjuster in values(path.adjusters)
					update_shift_width!(adjuster)
				end
			end
		end

		#generates observable for each cycle of "observable_skip"
		if mod(step, observable_skip) == 0 && step > equilibrium_skip
			#observables_counter += 1

			for observable in observables
				for estimator in estimators
					append!(output_observables[string(observable)][string(Symbol(estimator))], observable(path, potential, estimator))

				
				#output_observables[string(observable)][string(Symbol(estimator))] = (((observables_counter-1) * output_observables[string(observable)][string(Symbol(estimator))]) + observable(path, potential, estimator))/observables_counter  #rolling weighted mean of observable
				
				end
			end

		end
	end

	system_stats["acceptance_ratio"] = Dict()

	for mover in movers[1]
		system_stats["acceptance_ratio"][string(Symbol(mover))] = system_stats["acceptance_array"][string(Symbol(mover))] / system_stats["attempted_array"][string(Symbol(mover))] 
	end
	
	
	return [system_stats["acceptance_ratio"], output_observables]

end

function primitive_action(path::Path, bead::Int, particle::Int, potentials::Array{Potential})
    sum(kinetic_action(path, bead, particle), potential_action(path, bead, particle, potentials))
end

function PIMC(n_steps::Int, path::Path, movers, observables, potentials::Union{Potential, Array{Potential}}; n_particles = 1, visual = true)
	
	observable_skip = 0.001 * n_steps
	equilibrium_skip = 0.0 * n_steps
	
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

