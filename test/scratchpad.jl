using Revise # Tim Holy's Revise.jl; will auto-reload changed code.

ENV["EDITOR"] = "vim" # Natch.
push!(LOAD_PATH,"../src/") # load module from local directory

using HALCYON

V(x)=10*x^4 - x^2
HALCYON.sampledV(V)

