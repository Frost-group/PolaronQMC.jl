
using Plots
using Statistics


#Initialising


begin
#path
    mutable struct Path
        positions :: Array{Float64,1}
        n_beads :: Int
        τ :: Float64
        mass :: Float64
        n_dimensions :: Int
        λ :: Float64
        n_particles :: Int


        function Path(n_beads::Int, mass::Float64, n_dimensions::Int, n_particles::Int)
            τ = β/n_beads
            positions = rand(n_beads)*startrange*rand([-1,1])
            λ = ħ/2*mass 
            new(positions, n_beads, τ, mass, n_dimensions, λ, n_particles)
        end
    end
end


begin #Potentials
    begin #Potentials definition
        abstract type Potential end
        abstract type OneBodyPotential <: Potential end
        abstract type TwoBodyPotential <: Potential end

        struct HarmonicPotential <: OneBodyPotential
            ω :: Float64
            function HarmonicPotential(ω)
                new(ω)
            end
        end  
    end 
    begin #Potential functions
        function one_body_potential(potential::HarmonicPotential, path, pos)
            output = 0.5*path.mass*ω^2*pos^2 
            #output *= path.τ
            return output
        end
    end
end

begin #actions
    function kinetic_action(path::Path, pos1::Float64, pos2::Float64) #bead refers to index of current bead
        preterm = path.n_dimensions/2*path.n_particles*log(4*pi*path.λ*path.τ)
        return preterm + ((pos1-pos2)^2)/(4*path.λ*path.τ)
    end
    function action(path::Path, pos1::Float64, pos2::Float64, potential::Potential)
        return kinetic_action(path,pos1,pos2)+one_body_potential(potential,path,pos1)

    end
end

begin #move function
    function move(path::Path,bead::Int,h::Float64, potential::Potential)
        u = rand((-h,h))*rand()
        proposedpos = path.positions[bead] + u #proposed position to move to
        if bead < path.n_beads
            old_action = action(path, path.positions[bead], path.positions[bead+1], potential)
            new_action = action(path, proposedpos, path.positions[bead+1], potential)
        elseif bead == path.n_beads
            old_action = action(path, path.positions[bead], path.positions[1], potential)
            new_action = action(path, proposedpos, path.positions[1], potential)
        end

        δS = new_action-old_action
        prob = minimum([1,exp(-δS)])
        if rand()<prob
            path.positions[bead] = proposedpos
        end
    end
end


function Bisect(path,clip_length)
    sampled_beads = [Int(floor(clip_length/2))]
    to_cut = [Int(floor(clip_length/2))]
    for x in to_cut
        
        
    end
end


function Multi_Bisect!(path::Path, potential::Potential)
    max_level = Int(floor(rand(1:path.n_beads)/log(2)))
    for level in 1:max_level
        clip_length = 2^level - 1
    end
        
end


begin #PIMC
    function PIMC(path::Path, h::Float64, potential::Potential, niter::Int)
        for i in 1:niter
            site = rand(1:path.n_beads)
            move(path,site,h,potential)
        end
        
    end

end

#---------- Methods of sampling -------------------






#----------------------------------- Getting observables -----------------------
begin #Estimators
    abstract type Estimator end

    #Thermodynamic Estimator
    struct ThermodynamicEstimator <: Estimator end

end


begin #energy

    #Thermodynamic Estimator
    function Kinetic_energy(path, estimator::ThermodynamicEstimator)
        KE = 0
        term_one = (path.n_dimensions*path.n_particles) / (2*path.τ) #first term
        link_factor = 1/(4*path.λ*path.τ^2)
        #link_factor = 1
        for i in 1:path.n_beads
            if i>1
                link_term = ((path.positions[i]-path.positions[i-1])^2)*link_factor
            elseif i==1
                link_term = ((path.positions[i]-path.positions[path.n_beads])^2)*link_factor
            end

            KE += (term_one - link_term)
            
        end
        println("KE= ",KE)
        return KE
    end

        

    function Potential_energy(path::Path, potential::HarmonicPotential)
        PE = 0
        for i in 1:path.n_beads
            PE += 0.5*path.mass*ω^2*path.positions[i]^2
        end
        println("PE= ",PE)
        return PE
    end

    function Energy(path::Path, estimator::Estimator, potential::Potential)
        return (Kinetic_energy(path, estimator) + Potential_energy(path,potential))/path.n_beads
    end
end

begin #ground state wave function
    function wave_function(path::Path)
    end
end





#--------------------- comparison to actual
begin
    function analytic_energy(path::Path,potential::Potential,Max_n::Int)
        energy = 0
        partition_function = 0
        #partition function
        for i in 1:Max_n
            one_state_energy = ħ*ω*(i+0.5)
            partition_function += exp(-1*one_state_energy*β)
        end
        
        #av_energy
        for i in 1:Max_n
            one_state_energy = ħ*ω*(i+0.5)
            energy += one_state_energy * (exp(-one_state_energy*β)/partition_function)
                
        end
        return energy
    end
end








#Starting program (Init params)
begin
    startrange = 1
    ħ = 1
    k = 1

    n_beads = 100
    mass = 1.0
    ω = 1.0
    h = 3.0

    n_dimensions = 1
    n_particles = 1
    T = 1.0
    β = 1/(k*T)


    niter = 200000
    potential = HarmonicPotential(ω)
    estimator = ThermodynamicEstimator()



    path = Path(n_beads ,mass ,n_dimensions, n_particles)
    PIMC(path,h,potential,niter)

    average_pos = mean(path.positions)
    println("Average position is ",average_pos)



    println("System energy is ",Energy(path,estimator,potential))
    println("Analytic energy is ",analytic_energy(path,potential,10000))
#=
    plot(path.positions,1:path.n_beads)
    scatter!(path.positions,1:path.n_beads)
    xlabel!("Position")
    ylabel!("τ")
=#
    histogram(path.positions,bins=20)



end




begin
a = [1,2,3,4,5,6]
end



