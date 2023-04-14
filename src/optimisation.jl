# optimisation.jl


function updateAdjuster(adjuster::Union{SingleAdjuster, DisplaceAdjuster}, potential::Potential)
    if (adjuster.attempt_counter != 0)
        adjuster.acceptance_rate = adjuster.success_counter / adjuster.attempt_counter

        if adjuster.acceptance_rate < 0.01
            adjuster.value *= 0.99

        elseif adjuster.acceptance_rate > 0.99
            adjuster.value *=  1.01

        else
            adjuster.value *= (adjuster.acceptance_rate / 0.5) #originally it's 0.5
        end

        adjuster.attempt_counter = 0
        adjuster.success_counter = 0
    end
end


function updateAdjuster(adjuster::Union{BisectAdjuster}, potential::HarmonicPotential)

    adjuster.acceptance_rate = adjuster.success_counter / adjuster.attempt_counter

    if adjuster.acceptance_rate < 0.01
        adjuster.value *= 0.9

    elseif adjuster.acceptance_rate > 0.99
        adjuster.value *=  1.1

    else
        adjuster.value *= (adjuster.acceptance_rate / 0.5)
    end

    adjuster.attempt_counter = 0
    adjuster.success_counter = 0
end


function updateAdjuster(adjuster::Union{BisectAdjuster}, potential::FrohlichPotential)

    adjuster.acceptance_rate = adjuster.success_counter / adjuster.attempt_counter

    adjuster.attempt_counter = 0
    adjuster.success_counter = 0
end



function thermalised_start!(path::Path, potential::Potential; n_steps::Int = 2000, movers::Array = [[Bisect!],[1.0]], threads::Bool = true, verbose::Bool = true)
    st_regime = Primitive_Regime()
    st_observables = []
    st_estimators = []
    PIMC(n_steps, n_steps * 2 , n_steps * 2, path, movers, st_observables, st_estimators, potential, st_regime, adjust=true, threads=threads)
    if verbose
        println("Thermalisation complete")
    end
end






