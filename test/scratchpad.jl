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
ρ=zeros(evecs[:,1]) # size of discretised density
for i in 1:n
    ρ+=evecs[:,i].^2 
    plot!(evecs[:,i].^2)
end

plot!(ρ)

gui() #open dem plots
