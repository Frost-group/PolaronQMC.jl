#comparison.jl

#Real values of observables to compare against


function analytic_energy_harmonic(potential::HarmonicPotential, β::Float64, ħ::Float64)
    return potential.ω*ħ/2 + potential.ω*ħ/(exp(ħ*potential.ω*β)-1)
end


function average_position(path,particle)
    for i in 1:path.n_beads
           
