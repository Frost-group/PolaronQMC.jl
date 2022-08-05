



# function used to update shift width

function update_shift_width!(adjuster::Union{Single_Adjuster, Displace_Adjuster})
    if adjuster.adjust_counter >= 4
        #println("adjusted +") 

        adjuster.shift_width *= 1.2
        adjuster.adjust_counter = 0
    elseif adjuster.adjust_counter <= -4

        adjuster.shift_width /= 1.2
        adjuster.adjust_counter = 0
    end
end

function update_shift_width!(adjuster::Bisect_Adjuster)
    for level in 1:adjuster.max_level 
        if adjuster.adjust_counter_array[string(level)] >= 5
            #println("adjusted +")
            adjuster.shift_width_array[string(level)] *= 1.2
            adjuster.adjust_counter_array[string(level)] = 0


        elseif adjuster.adjust_counter_array[string(level)] <= -5
            #println("adjusted -")
            adjuster.shift_width_array[string(level)] /= 1.2
            adjuster.adjust_counter_array[string(level)] = 0
        end
    end
end








