



"""
update_shift_width!(adjuster::Union{Single_Adjuster, Displace_Adjuster})

Updates the "shift width" parameter used in Adjusters within Movers.

# Arguments
- `adjuster`: Type of adjuster specific for the type of move being attempted.


"""

function update_shift_width!(adjuster::Union{Single_Adjuster, Displace_Adjuster})

    adjuster.acceptance_rate = adjuster.success_counter / adjuster.attempt_counter

    if adjuster.acceptance_rate < 0.01
        adjuster.shift_width *= 0.99

    elseif adjuster.acceptance_rate > 0.99
        adjuster.shift_width *=  1.01

    else
        adjuster.shift_width *= (adjuster.acceptance_rate / 0.5)
    end

    adjuster.attempt_counter = 0
    adjuster.success_counter = 0
end


function update_shift_width!(adjuster::Union{Bisect_Adjuster})

    adjuster.acceptance_rate = adjuster.success_counter / adjuster.attempt_counter

    #if adjuster.acceptance_rate < 0.01
        #adjuster.value *= 0.9

    #elseif adjuster.acceptance_rate > 0.99
        #adjuster.value *=  1.1

    #else
        #adjuster.value *= (adjuster.acceptance_rate / 0.5)
    #end

    adjuster.attempt_counter = 0
    adjuster.success_counter = 0
end

#=
function update_shift_width!(n_beads, adjuster::Bisect_Adjuster)

    adjuster.acceptance_rate = adjuster.success_counter / adjuster.attempt_counter
    
    if adjuster.acceptance_rate < 0.6
        if adjuster.value > 1
            adjuster.value -= 1
        end

    elseif adjuster.acceptance_rate > 0.6
        if 2^(adjuster.value+1)+1 <= n_beads
            adjuster.value += 1
        end
    end

    adjuster.attempt_counter = 0
    adjuster.success_counter = 0
end
=#


#=
function update_shift_width!(adjuster::Union{Single_Adjuster, Displace_Adjuster})

    if adjuster.success_counter < 5
        adjuster.shift_width = adjuster.shift_width * 1.1

    elseif adjuster.attempt_counter - adjuster.success_counter < 5
        adjuster.shift_width = adjuster.shift_width * 0.9

    else
        adjuster.acceptance_rate = adjuster.success_counter / adjuster.attempt_counter
        adjuster.shift_width = adjuster.shift_width * (0.234 / adjuster.acceptance_rate)
        adjuster.attempt_counter = 0
        adjuster.success_counter = 0       
    end
end

function update_shift_width!(adjuster::Union{Single_Adjuster, Displace_Adjuster})
    if adjuster.adjust_counter >= 5
        #println("adjusted +") 

        #weight average for new shift width
        adjuster.shift_width = (adjuster.shift_width*6 + adjuster.shift_width+adjuster.adjust_unit*4) / 10
        adjuster.adjust_counter = 0
    elseif adjuster.adjust_counter <= -5

        adjuster.shift_width = (adjuster.shift_width*6 - adjuster.shift_width+adjuster.adjust_unit*4) / 10
        adjuster.adjust_counter = 0
    end
end

function update_shift_width!(adjuster::Bisect_Adjuster)
end

=#


"""
    thermalised_start!(path::Path, potential::Potential; n_steps::Int = 2000, movers::Array = [[Bisect!],[1.0]])

Performs PIMC simulation on the path without producing observables, thermalising the system before sampling of observables.


"""
function thermalised_start!(path::Path, potential::Potential; n_steps::Int = 2000, movers::Array = [[Bisect!],[1.0]], threads::Bool = true, verbose::Bool = true)
    st_regime = Primitive_Regime()
    st_observables = []
    st_estimators = []
    PIMC(n_steps, n_steps * 2 , n_steps * 2, path, movers, st_observables, st_estimators, potential, st_regime, adjust=true, threads=threads)
    if verbose
        println("Thermalisation complete")
    end
end






