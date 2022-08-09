



# function used to update shift width

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


function thermalised_start!(path::Path, n_steps::Int, movers::Array, potential::Potential)
    st_regime = Primitive_Regime()
    st_observables = []
    st_estimators = []
    pimc = PIMC(n_steps, n_steps * 2 , n_steps * 2, path, movers, st_observables, st_estimators, potential, st_regime, adjust=true)
end






