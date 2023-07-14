#errors.jl
using StatsBase
using Statistics


function jackknife(observables_array)
    """
    Jackknife error analysis based on "User's Guide to Monte Carlo Methods for evaluating path integrals"
    (Dividing data into blocks of bins)
    """

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

    # Getting jackknife estimators
    jk_estimators = []

    observables_sum = sum(observables_array)
    for i in 1:n_bins
        jk_estimator = observables_sum - (bin_width * block_estimators[i])
        jk_estimator /= jk_binwidth

        append!(jk_estimators,jk_estimator)
    end

    # Getting jackknife variance
    variance_jk = 0.0
    for k in 1:n_bins
        variance_jk += (jk_estimators[k] - mean(observables_array))^2
    end
    variance_jk *= (n_bins-1)/n_bins

    return [variance_bins,variance_jk]
end

function autoCorrelation(observable_arr)
    """
    Autocorrelation function to determine if the data are taken in time intervals beyond the autocorrelation time 
    
    Parameter:
        observable_arr
    """

    # Length of autocorrelation is equal to the total number of observables entries minus 1
    autoCorrelation = zeros(length(observable_arr)-1)

    # Expected value and variance
    O_avg = mean(observable_arr)
    len = length(observable_arr)
    Var = var(observable_arr)

    # Calculating autocorrelation coefficients
    for k in 1:length(observable_arr)-1

        for i in 1:length(observable_arr)
            k2 = k + i
            if k2 > length(observable_arr)
                break
            end
            E1 = observable_arr[i]
            E2 = observable_arr[k2]
            autoCorrelation[k] += (E1 - O_avg) * (E2 - O_avg) / (len-k)
        end
        autoCorrelation[k] /= Var 
    end

    # Return an array of autocorrelation coefficients
    # The faster the decay, the quicker the de-correlation process
    return autoCorrelation
end

function autoCorrelationTime(Ck)
    """
    return the autocorrelation time based on autocorrelation coefficients
    Method based on Ceperly "interacting electrons"
    """
    return 1 + 2*sum(Ck)
end



    
