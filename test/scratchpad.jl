push!(LOAD_PATH,"../src/") # load module from local directory

using HALCYON

V(x)=10*x^4 - x^2
HALCYON.sampledV(V)

