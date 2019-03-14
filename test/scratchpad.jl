using Revise # Tim Holy's Revise.jl; will auto-reload changed code.

ENV["EDITOR"] = "vim" # Natch.
push!(LOAD_PATH,"../src/") # load module from local directory

using HALCYON

println("Loading potential...")
V=HALCYON.DoubleWell # our 1D potential.
println("Sampled values...")
HALCYON.sampledV(V)

println("Solving Time Indie. Schr. Eqn.")
evals,evecs=HALCYON.TISE(V)

using Plots
plot(HALCYON.sampledV(V))

# calc density for a set of (non-interacting) states
n=4
ρ=similar(evecs[:,1]) # size of discretised density
fill!(ρ,0)

for i in 1:n
    ρ+=evecs[:,i].^2 
    plot!(evecs[:,i].^2)
end

plot!(ρ)

gui() #open dem plots

using GaussianProcesses

# https://github.com/STOR-i/GaussianProcesses.jl/blob/master/notebooks/Regression.ipynb
# Training data

#Select mean and covariance function
mZero = MeanZero()                   #Zero mean function
kern = SE(0.0,0.0)                   #Sqaured exponential kernel (note that hyperparameters are on the log scale)
logObsNoise = -1.0                        # log standard deviation of observation noise (this is optional)

# OK; pretty coarse idea! Just fit density to bare potential. :^) 
gp = GP(HALCYON.sampledV(V),ρ,mZero,kern,logObsNoise)       #Fit the GP

using Plots
plot(gp)
# Ah, magic!
