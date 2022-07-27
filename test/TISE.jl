using Revise # Tim Holy's Revise.jl; will auto-reload changed code.

ENV["EDITOR"] = "vim" # Natch.
push!(LOAD_PATH,"../src/") # load module from local directory

using Halcyon

println("Loading potential...")
V=Halcyon.DoubleWell # our 1D potential.
println("Sampled values...")
Halcyon.sampledV(V)

println("Solving Time Indie. Schr. Eqn.")
evals,evecs=Halcyon.TISE(V)

using Plots
plot(Halcyon.sampledV(V))

# calc density for a set of (non-interacting) states
n=4
ρ=similar(evecs[:,1]) # size of discretised density
fill!(ρ,0)

for i in 1:n
    global ρ+=evecs[:,i].^2 
    plot!(evecs[:,i].^2)
end

plot!(ρ)

gui() #open dem plots

