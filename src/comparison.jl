#comparison.jl


#Real values of observables to compare against


function analyticEnergyHarmonic(ω::Float64, β::Float64, ħ::Float64, n_dimensions::Int)
    """
        Analytic average energy <E> for a harmonic oscillator
    """
    return n_dimensions / 2 * ħ * ω * (1 .+ exp.(-ħ * ω * β)) ./ (1 .- exp.(-ħ * ω * β))
end


function selectiveMean(observable_array::Array, limit::Union{Int,Float64})
    """
    mean that filters out extreme values
    """
    counter = 0
    mean_value = 0.0
    for observable_value in observable_array
        if observable_value > abs(mean_value + 1) * -limit &&
           observable_value < abs(mean_value + 1) * limit
            counter += 1

            mean_value = ((counter - 1) * mean_value + 1 * observable_value) / counter

        end
    end
    return mean_value
end
