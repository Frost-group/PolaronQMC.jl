#errors.jl
using StatsBase
using Statistics

#methods for finding the errors on the estimator

function jackknife(observables_array)
    bin_width = Int(floor(0.1*length(observables_array))) #block binwidth
    jk_binwidth = length(observables_array)-bin_width #jack knife binwidth, for complementary bins
    n_bins = cld(length(observables_array),bin_width)


    #getting block estimators

    block_bins = collect(Iterators.partition(observables_array,bin_width)) 
    block_estimators = [sum(bin)/bin_width for bin in block_bins]

    bins = []
    for i in 1:n_bins
        append!(bins,[sample(observables_array,jk_binwidth,replace=false)])
    end

    #getting jack knife estimators
    jk_estimators = []
    for i in 1:length(bins)
        jk_estimator = 0
        for observable in bins[i]
            jk_estimator += observable - bin_width * block_estimators[i]
        end
        jk_estimator /= jk_binwidth
        append!(jk_estimators,jk_estimator)
    end

    #getting jack knife variance
    jk_variance = 0.0
        for k in 1:n_bins
            jk_variance += (jk_estimators[k] - mean(observables_array))^2
        end
        jk_variance *= (n_bins-1)/n_bins
        return jk_variance

end


    
