begin
using Revise
using PolaronQMC
using BenchmarkTools



#variables used
    T = 3.0
    β = 1/T
    alpha_range = 1.0
    thermalisation_steps = 100
    steps_base = 800
    fixed_τ = 0.001
    adjusted_beads = Int(floor(1/(fixed_τ*T)))

    path = Path(adjusted_beads, 1, 3, fixed_τ)
    particle = rand(1:path.n_particles)
    bead = rand(1:path.n_beads)

    α = 1.0

    frohlich_potential = FrohlichPotential(α,1.0,1.0)
    harmonic_potential = HarmonicPotential(1.0)


    frohlichcache = PotentialCache(path, frohlich_potential)
    harmoniccache = PotentialCache(path, harmonic_potential)

    println("Benchmarking PIMC moves...")
    println("Beads: $adjusted_beads ")

    
end



begin
n_steps = 5000

println("PIMC using updated Bisect")
pimc_bm = @btime PIMC(n_steps, n_steps, n_steps, path, [[Bisect!],[1.0]], [], [], FrohlichPotential(α, 1.0, 1.0), Primitive_Regime(), threads = false)
println("\n")

println("Old PIMC")
pimc_bm = @btime PIMC(n_steps, n_steps, n_steps, path, [[BisectL!],[1.0]], [], [], FrohlichPotential(α, 1.0, 1.0), Primitive_Regime())
println("\n")


end
begin


#Single


println("Single!() Frohlich")
single_bm = @btime Single!(path, particle, frohlich_potential, frohlichcache, Primitive_Regime(), path.adjusters["Single!"])
println("\n")

println("SingleL!() Frohlich")
singleL_bm = @btime SingleL!(path, particle, frohlich_potential, frohlichcache, Primitive_Regime(), path.adjusters["Single!"])
println("\n")



#Bisect
println("Bisect!() Harmonic")
bisect_bm = @btime Bisect!(path, particle, harmonic_potential, harmoniccache, Primitive_Regime(), path.adjusters["Bisect!"])
println("\n")

println("Bisect!() Frohlich")
bisect_bm = @btime Bisect!(path, particle, frohlich_potential, frohlichcache, Primitive_Regime(), path.adjusters["Bisect!"])
println("\n")


println("BisectL!() Frohlich")
bisectL_bm = @btime BisectL!(path, particle, frohlich_potential, frohlichcache, Primitive_Regime(), path.adjusters["Bisect!"])
println("\n")


end
begin
#benchmarking smaller sections


println("Updated Frohlich one_body_potential")
obp_bm = @btime one_body_potential(frohlich_potential, frohlichcache, path, bead, particle)
println("\n")

println("Previous Frohlich one_body_potential")
obp_bm = @btime one_body_potentialL(frohlich_potential, frohlichcache, path, bead, particle)
println("\n")


end
