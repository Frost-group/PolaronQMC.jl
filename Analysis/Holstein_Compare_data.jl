# Holstein Polaron Perturbation theory and DMC data from literature
# Data extracted from https://amslaurea.unibo.it/22104/1/ragni-diagrammatic_monte_carlo_study_of_the_holstein_polaron.pdf

function PT_Holstein_Strong(alpha, dim)
    # t = 1.0, ω = 1.0
    return -2*exp.(-alpha*dim)*dim .- alpha * dim .- 1 ./ alpha
end

function PT_Holstein_Weak(alpha, dim)
    # t = 1.0, ω = 1.0
    return -2 * dim .- alpha .* 2.8099/(2π)
end

DMC_α_1D = [0.0, 0.25, 0.5, 0.75, 1.0, 1.3, 1.57, 1.84, 2.1, 2.36
            , 2.61, 2.88, 3.15, 3.42, 3.67, 3.94, 4.21, 4.47, 5.0]
DMC_α_E_1D = [-2.0, -2.12, -2.25, -2.37, -2.51, -2.63, -2.78, -2.91, -3.07, -3.2, -3.4, 
            -3.54, -3.73, -3.93, -4.12, -4.33, -4.56, -4.78, -5.25]

DMC_α_2D = [0.0, 0.25, 0.51, 0.78, 1.05, 1.3, 1.57, 1.84, 2.1, 2.36
            , 2.63, 2.9, 3.15, 3.42, 3.68, 3.95, 4.21, 4.48, 4.74, 5.0]
DMC_α_E_2D = [-4.00, -4.14, -4.28, -4.44, -4.58, -4.75, -4.92, -5.11, -5.31, -5.58, -5.88,
                -6.25, -6.71, -7.17, -7.69, -8.15, -8.66, -9.16, -9.67, -10.19]

DMC_α_3D = [0.0, 0.26, 0.53, 0.8, 1.07, 1.33, 1.6, 1.86, 2.13, 2.4, 2.66, 2.93, 3.2, 3.46,
            3.73, 4.0]
DMC_α_E_3D = [-6.00, -6.13, -6.27, -6.42, -6.57, -6.74, -6.91, -7.11, -7.36, -7.76, -8.44,
            -9.21, -9.98, -10.73, -11.54, -12.30]
