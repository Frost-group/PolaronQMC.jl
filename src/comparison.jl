#comparison.jl


#Real values of observables to compare against


function analytic_energy_harmonic(potential::HarmonicPotential, β::Float64, ħ::Float64, n_dimensions::Int)
    return n_dimensions/2*ħ*potential.ω * (1 .+ exp.(-ħ*potential.ω*β))./(1 .- exp.(-ħ*potential.ω*β))
    #return potential.ω*ħ/2 + potential.ω*ħ/(exp(ħ*potential.ω*β)-1)
end


function analytic_energy_harmonic2(potential::HarmonicPotential, β::Float64, ħ::Float64)
    return potential.ω*ħ/2 * (1+exp(-potential.ω*ħ*β))/(1-exp(-ħ*potential.ω*β))
end


function selective_mean(observable_array::Array, limit::Union{Int,Float64})
    """
    mean that filters out extreme values
    """
    counter = 0
    mean_value = 0.0
    for observable_value in observable_array
        if observable_value > abs(mean_value + 1)*-limit && observable_value < abs(mean_value + 1)*limit
            counter += 1
   
            mean_value = ((counter-1) * mean_value + 1 * observable_value) / counter

        end
    end
    return mean_value
end





           
