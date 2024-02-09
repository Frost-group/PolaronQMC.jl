using BenchmarkTools
using LinearAlgebra
using TimerOutputs
using StaticArrays
using PolaronQMC

#const tmr1 = TimerOutput();
#reset_timer!(tmr1)

begin
    n_beads = 100
    path = Path(n_beads, 1, 3, 1.0 / (0.1 * n_beads), m = 1.0)
    mover = SingleMover(path)
    potential = FrohlichPotential(4.0, 1.0, 1.0) # ω is phonon frequency
    #potential = HarmonicPotential(1.0)
    Action_arr = [0.0, 0.0]
    shift = zeros(3)
    regime = PrimitiveRegime()
    x = zeros(3)
    particle = 1
    path1 = deepcopy(path)
    path2 = deepcopy(path)
end
#=
@time begin
    particle = 1
    @timeit tmr1 "Moving" for sweep in 1:10000
        moveBead!(mover, path, particle, potential, regime, Action_arr, shift, well_size = 4.0) # Moving beads a total of n_sweep times
    end
end
=#
begin
    fac = 1.0 * path.τ * path.n_beads * 1.0
    #c = Array{Float64}(undef, n_beads, n_beads)
    c = zeros(n_beads, n_beads)
    for i = 1:n_beads
        for j = i+1:n_beads
            c[j, i] = cosh(fac * ((j - i) / n_beads - 0.5))
        end
    end
end

@time begin
    #@timeit tmr1 "Potential" 
    for _ = 1:200000
        bead = rand(1:n_beads)
        2 * path.τ * oneBodyPotential(potential, path1, bead, particle, x, c)
    end
end

@time begin
    #@timeit tmr1 "Potential"
    d = zeros(n_beads, n_beads)
    temp = zeros(3)
    for i = 1:n_beads
        for j = i+1:n_beads
            d[j, i] = norm([
                path2.beads[j, particle, k] - path2.beads[i, particle, k] for
                k = 1:path2.n_dimensions
            ])
        end
    end
    for _ = 1:200000
        bead = rand(1:n_beads)
        2 * path.τ * oneBodyPotential2(potential, path2, bead, particle, c, x, d)
    end
end
#println(tmr1)
#reset_timer!(tmr1)
#=
@time begin 
    #d = Array{Float64}(undef, 5, 5)
    d = zeros(5, 5)
    for i in 1:5, j in i+1:5
        d[j, i] = cosh((i-j)/100-0.5)
    end
    #println(d)
end
=#
#=
@btime begin 
    #d = Array{Float64}(undef, 5, 5)
    d = zeros(5, 5)
    for i in 1:5, j in i:5
        d[i, j] = cosh((i-j)/100-0.5)
    end
   #println(d)
end
=#
