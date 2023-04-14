#errors.jl
using StatsBase
using Statistics


function jackknife(observables_array)

    # Block binwidth
    bin_width = Int(floor(0.1*length(observables_array))) 
    jk_binwidth = length(observables_array)-bin_width 
    # Jack knife binwidth, for complementary bins
    n_bins = cld(length(observables_array),bin_width)

    # Getting block estimators
    block_bins = collect(Iterators.partition(observables_array,bin_width)) 
    block_estimators = [sum(bin)/bin_width for bin in block_bins]

    # Getting bin_variance
    variance_bins = 0.0
    for k in 1:n_bins
        variance_bins += (block_estimators[k] - mean(observables_array))^2
    end
    variance_bins /= 1/(n_bins*(n_bins-1))

    # Getting jack knife estimators
    jk_estimators = []

    observables_sum = sum(observables_array)
    for i in 1:n_bins
        jk_estimator = observables_sum - (bin_width * block_estimators[i])
        jk_estimator /= jk_binwidth

        append!(jk_estimators,jk_estimator)
    end

    # Getting jack knife variance
    variance_jk = 0.0
    for k in 1:n_bins
        variance_jk += (jk_estimators[k] - mean(observables_array))^2
    end
    variance_jk *= (n_bins-1)/n_bins

    return [variance_bins,variance_jk]
end

function autoCorrelation(observable_arr, observable_skip)
    len = length(observable_arr)
    """
    observable_arr: Usually energy array
    """
    autoCorrelation = zeros(length(observable_arr)-1)
    O_avg = mean(observable_arr)

    Var = var(observable_arr) * (len-1)

    for i in 1:length(observable_arr)-1

        for k in 1:length(observable_arr)
            k2 = i + k
            if k2 > length(observable_arr)
                break
            end
            E1 = observable_arr[k]
            E2 = observable_arr[k2]
            autoCorrelation[i] += (E1-O_avg) * (E2 - O_avg)
        end
        autoCorrelation[i] /= Var
    end
    return autoCorrelation
end

function autoCorrelationTime(Ck)
    return 1 + 2*sum(Ck)
end



    
