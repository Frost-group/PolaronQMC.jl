

function update_shift_width!(adjuster::Union{Single_Adjuster, Displace_Adjuster}, potential::HarmonicPotential)

    adjuster.acceptance_rate = adjuster.success_counter / adjuster.attempt_counter

    if adjuster.acceptance_rate < 0.01
        adjuster.value *= 0.99

    elseif adjuster.acceptance_rate > 0.99
        adjuster.value *=  1.01

    else
        adjuster.value *= (adjuster.acceptance_rate / 0.5)
    end

    adjuster.attempt_counter = 0
    adjuster.success_counter = 0
end


function update_shift_width!(adjuster::Union{Bisect_Adjuster}, potential::HarmonicPotential)

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


function update_shift_width!(adjuster::Union{Bisect_Adjuster}, potential::FrohlichPotential)

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






