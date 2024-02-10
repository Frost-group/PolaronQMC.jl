# PolaronQMCVisualisation.jl
using Plots
using LaTeXStrings


function draw_beads_3d(beads, xlims, ylims, zlims, n_particles, frame, potential, mover, T)

    default(
        fontfamily = "Computer Modern",
        titlefont = (16, "Computer Modern"),
        legendfontsize = 12,
        guidefont = (18, "Computer Modern"),
        tickfont = (12, "Computer Modern"),
        linewidth = 2,
        framestyle = :box,
        label = nothing,
        grid = true,
    )

    p = plot(xlabel = L"x", ylabel = L"y", zlabel = L"z")

    for particle = 1:n_particles

        x = beads[:, particle, 1]
        y = beads[:, particle, 2]
        z = beads[:, particle, 3]

        x = Array(reshape(x, length(x)))
        y = Array(reshape(y, length(y)))
        z = Array(reshape(z, length(z)))

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
        title!(
            "PIMC Animation [Frame=$(frame)] \n $(potential) \n $(mover) \n T=$(T)",
            title_position = :left,
        )

    end
    return p
end


function draw_beads_2d(beads, xlims, ylims, n_particles, frame, potential, mover, T)

    default(
        fontfamily = "Computer Modern",
        titlefont = (16, "Computer Modern"),
        legendfontsize = 12,
        guidefont = (18, "Computer Modern"),
        tickfont = (12, "Computer Modern"),
        linewidth = 2,
        framestyle = :box,
        label = nothing,
        grid = true,
    )

    p = plot(xlabel = L"x", ylabel = L"y")

    for particle = 1:n_particles

        x = beads[:, particle, 1]
        y = beads[:, particle, 2]

        x = Array(reshape(x, length(x)))
        y = Array(reshape(y, length(y)))

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
        title!(
            "PIMC Animation [Frame=$(frame)] \n $(potential) \n $(mover) \n T=$(T)",
            title_position = :left,
        )

    end
    return p
end


function animate_PIMC(data, n_particles, n_dimensions, potential, mover, T)

    if n_dimensions == 3
        xlims = [
            minimum([minimum(x[:, :, 1]) for x in data[3]]),
            maximum([maximum(x[:, :, 1]) for x in data[3]]),
        ]
        ylims = [
            minimum([minimum(x[:, :, 2]) for x in data[3]]),
            maximum([maximum(x[:, :, 2]) for x in data[3]]),
        ]
        zlims = [
            minimum([minimum(x[:, :, 3]) for x in data[3]]),
            maximum([maximum(x[:, :, 3]) for x in data[3]]),
        ]

        frame = 0
        animation = Plots.@animate for p in data[3]
            draw_beads_3d(p, xlims, ylims, zlims, n_particles, frame, potential, mover, T)
            frame += 1
        end

    elseif n_dimensions == 2
        xlims = [
            minimum([minimum(x[:, :, 1]) for x in data[3]]),
            maximum([maximum(x[:, :, 1]) for x in data[3]]),
        ]
        ylims = [
            minimum([minimum(x[:, :, 2]) for x in data[3]]),
            maximum([maximum(x[:, :, 2]) for x in data[3]]),
        ]

        frame = 0
        animation = Plots.@animate for p in data[3]
            draw_beads_2d(p, xlims, ylims, n_particles, frame, potential, mover, T)
            frame += 1
        end
    end

    return animation
end

function animate_PIMC(data, n_dimensions; particle = 1, dim = 1)
    positions = data["Position:p$(particle)d$(dim)"]
    positions_flatten = collect(Iterators.flatten(positions))
    animation = Plots.@animate for i = 1:length(positions)
        scatter(
            positions[i],
            1:length(positions[1]),
            xlimit = [minimum(positions_flatten), maximum(positions_flatten)],
            xlabel = "x",
            ylabel = "τ",
        )
    end

    return animation
end
