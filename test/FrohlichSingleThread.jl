
pot = "Frohlich";
T = 0.5; version=rand(1:10000)
data, energy, error, comparison_energy, equilibrium_skip, observables_skip, n_beads, path = 
    generalPIMC(T, 1.0, 2.0, 3, 100000, version=1, pot=pot)

println("-----Simulation Ended-----")


