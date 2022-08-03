

#Adjuster for the Single! move algorithm
mutable struct Single_Adjuster <: Adjuster
    adjust_counter :: Int
    shift_width :: Float64
    adjust_unit :: Float64 #how much shift width is adjusted by each time
    function Single_Adjuster(path::Path)
        shift_width = sqrt(4 * path.λ * path.τ)
        adjust_unit = shift_width
        new(0,shift_width, adjust_unit)
    end
end


#Adjuster for the Displace! move algorithm
mutable struct Displace_Adjuster <: Adjuster
    adjust_counter :: Int
    shift_width :: Float64 
    adjust_unit :: Float64 #how much shift width is adjusted by each time
    function Displace_Adjuster(path::Path)
        shift_width = sqrt(4 * path.λ * path.τ)*10
        adjust_unit = 0.5*shift_width
        new(0,shift_width, adjust_unit)
    end
end




#Adjuster for the Bisect! move alogrithm
mutable struct Bisect_Adjuster <: Adjuster
    adjust_counter_array :: Dict
    shift_width_array :: Dict
    max_level :: Int


    function Bisect_Adjuster(path::Path)
        adjust_counter_array = Dict()
        shift_width_array = Dict()
        
        #max_level = Int(floor(log(rand(1:path.n_beads)) / log(2)))
        max_level = 4

        for level in 0:max_level
            shift_width_array[string(level)] = sqrt(2^(level) * path.λ * path.τ )
            adjust_counter_array[string(level)] = 0
        end
        
        new(adjust_counter_array,shift_width_array, max_level)
    end
end


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








