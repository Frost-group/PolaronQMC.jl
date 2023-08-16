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
    fixed_τ = 0.0025
    adjusted_beads = Int(floor(1/(fixed_τ*T)))

    path = Path(adjusted_beads, 1, 3, fixed_τ)

    α = 1.0


    println("Benchmarking PIMC moves...")
    println("Beads: $adjusted_beads ")

    particle = rand(1:path.n_particles)
    bead = rand(1:path.n_beads)


end


begin

println("PIMC")
pimcx_bm = @btime PIMC(100, 100, 100, path, [[Bisect!],[1.0]], [], [], FrohlichPotential(α, 1.0, 1.0), Primitive_Regime())
println("\n")

println("PIMCX")
pimcx_bm = @btime PIMCX(100, 100, 100, path, [[Bisect!],[1.0]], [], [], FrohlichPotential(α, 1.0, 1.0), Primitive_Regime())
println("\n")

end
begin


#Single
println("Single!()")
single_bm = @btime Single!(path, particle, FrohlichPotential(α, 1.0, 1.0), Primitive_Regime(), path.adjusters["Single!"])
println("\n")

#Bisect
println("Bisect!() Harmonic")
bisect_bm = @btime Bisect!(path, particle, HarmonicPotential(1.0), Primitive_Regime(), path.adjusters["Bisect!"])


println("Bisect!() Frohlich")
bisect_bm = @btime Bisect!(path, particle, FrohlichPotential(α, 1.0, 1.0), Primitive_Regime(), path.adjusters["Bisect!"])
println("\n")



#benchmarking smaller sections

println("Frohlich one_body_potential")
frohlich_bm = @btime one_body_potential(FrohlichPotential(α, 1.0, 1.0), path, bead, particle)
println("\n")

println("Deepcopy")
deepcopy_bm = @btime deepcopy(path.beads[:,particle, :])
println("\n")


end