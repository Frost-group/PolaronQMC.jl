# actions.jl



#Ways of calculating action

    #simple method


    function kinetic_action(path::Path, bead_one::Int, bead_two::Int, particle::Int, regime::Simple_Regime)
        kinetic_action = 0.5*path.m*(path.beads[bead_two, particle] - path.beads[bead_one, particle])^2
        return kinetic_action
    end


    function potential_action(path::Path, bead::Int, particle::Int, potential::OneBodyPotential, regime::Simple_Regime)
        return one_body_potential(potential, path, bead, particle)
    end


    function total_action(path::Path, bead_one::Int, bead_two::Int, particle::Int, potential::OneBodyPotential, regime::Simple_Regime)
        return kinectic_action(path, bead_one, bead_two, particle, regime) + potential_action(path, bead_one, particle, potential, regime)
    end









    #Primitive method (based off Ceperly paper)


    function kinetic_action(path::Path, bead_one::Int, bead_two::Int, particle::Int, regime::Primitive_Regime)
        kinetic_action = (path.beads[bead_two, particle] - path.beads[bead_one, particle])^2 / (4 * path.λ * path.τ)
        kinetic_action += 1.0 * path.n_particles / 2.0 * log(4π * path.λ * path.τ) #used to have n beads term, why?
        return kinetic_action
    end

    function potential_action(path::Path, bead::Int, particle::Int, potential::ConstantPotential, regime::Primitive_Regime)
        return path.τ * potential.V
    end

    function potential_action(path::Path, bead::Int, particle::Int, potential::OneBodyPotential, regime::Primitive_Regime)
        return path.τ * one_body_potential(potential, path, bead, particle)
    end

    function potential_action(path::Path, bead::Int, particle::Int, potential::TwoBodyPotential, regime::Primitive_Regime)
        potential_action = sum(two_body_potential(potential, path, bead, particle, other_particle) for other_particle in 1:path.n_particles if particle != other_particle)
        return path.τ * potential_action
    end


    function total_action(path::Path, bead_one::Int, bead_two::Int, particle::Int, potential::Potential, regime::Primitive_Regime)
        return kinetic_action(path, bead_one, bead_two, particle, regime) + potential_action(path, bead_two, particle, potential, regime)
    end

