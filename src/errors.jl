#errors.jl
using StatsBase
using Statistics

#methods for finding the errors on the estimator

function jackknife(observables_array)
    bin_width = Int(floor(0.1*length(observables_array))) #block binwidth
    # println("binwidth is: ", bin_width)
    jk_binwidth = length(observables_array)-bin_width #jack knife binwidth, for complementary bins
    n_bins = cld(length(observables_array),bin_width)

    #getting block estimators

    block_bins = collect(Iterators.partition(observables_array,bin_width)) 
    block_estimators = [sum(bin)/bin_width for bin in block_bins]

    #=
    bins = []
    for i in 1:n_bins
        append!(bins,[sample(observables_array,jk_binwidth,replace=false)])
    end
    =#

    #getting bin_variance

    variance_bins = 0.0
    for k in 1:n_bins
        variance_bins += (block_estimators[k] - mean(observables_array))^2
    end
    variance_bins /= 1/(n_bins*(n_bins-1))
    

    #getting jack knife estimators
    jk_estimators = []

    observables_sum = sum(observables_array)
    for i in 1:n_bins
        jk_estimator = observables_sum - (bin_width * block_estimators[i])
        jk_estimator /= jk_binwidth

        append!(jk_estimators,jk_estimator)
    end

    #getting jack knife variance
    variance_jk = 0.0
        for k in 1:n_bins
            variance_jk += (jk_estimators[k] - mean(observables_array))^2
        end
        variance_jk *= (n_bins-1)/n_bins
    return [variance_bins,variance_jk]
end



    
