



"""
    update_shift_width!(adjuster::Union{Single_Adjuster, Displace_Adjuster})

Updates the "shift width" parameter used in Adjusters within Movers.

# Arguments
- `adjuster`: Type of adjuster specific for the type of move being attempted.


"""

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
    for level in 1:adjuster.max_level 
        if adjuster.adjust_counter_array[string(level)] >= 7
            #println("adjusted +")

            #weight average for new shift width
            adjuster.shift_width_array[string(level)] = (adjuster.shift_width*6 + adjuster.shift_width+adjuster.adjust_unit*4) / 10
            adjuster.adjust_counter_array[string(level)] = 0


        elseif adjuster.adjust_counter_array[string(level)] <= -7
            #println("adjusted -")
            adjuster.shift_width_array[string(level)] = (adjuster.shift_width*6 - adjuster.shift_width+adjuster.adjust_unit*4) / 10
            adjuster.adjust_counter_array[string(level)] = 0
        end
    end
end



"""
    thermalised_start!(path::Path, potential::Potential; n_steps::Int = 2000, movers::Array = [[Bisect!],[1.0]])

Performs PIMC simulation on the path without producing observables, thermalising the system before sampling of observables.


"""
function thermalised_start!(path::Path, potential::Potential; n_steps::Int = 2000, movers::Array = [[Bisect!],[1.0]], threads::Bool = true)
    st_regime = Primitive_Regime()
    st_observables = []
    st_estimators = []
    pimc = PIMC(n_steps, n_steps * 2 , n_steps * 2, path, movers, st_observables, st_estimators, potential, st_regime, adjust=true, threads=threads)
    println("Thermalisation complete")
end






