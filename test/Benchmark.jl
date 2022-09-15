
using Revise
using PolaronQMC
using BenchmarkTools
println("Using ", Threads.nthreads(), " threads")

#variables used
T = 3.0
β = 1/T
alpha_range = 1.0
thermalisation_steps = 0
steps_base = 300
fixed_τ = 0.01
adjusted_beads = Int(floor(1/(fixed_τ*T)))
        
#benchmarking
total_steps = Int(sum(steps_base*alpha_range) + thermalisation_steps)



println("Benchmarking simulation...")
println("Running with $adjusted_beads beads")
println("total_steps = ", total_steps)

println(" ")
quickrun_1 = @benchmark quickrun_frohlich(T, alpha_range, fixed_τ, thermalisation_steps, steps_base, thermalised = false, verbose = false, threads=true)
println("Multi-threaded:")
println("Steps per second = ", mean(total_steps / (quickrun_1.times * 1e-9)))


println(" ")
quickrun_2 = @benchmark quickrun_frohlich(T, alpha_range, fixed_τ, thermalisation_steps, steps_base, thermalised = false, verbose = false, threads=false)
println("Single-threaded:")
println("Step per second = ", mean(total_steps / (quickrun_2.times * 1e-9)))




#show(io,MIME("text/plain"),b)







