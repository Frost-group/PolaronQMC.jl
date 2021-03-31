# PIMC.jl
# Initially location for all the path integral bits and bobs

# Very generic path
struct Path
    nbeads ::Int64
    nparticles ::Int64
    mass ::Float64
    spatialdimensions ::Int64

    R::Array{Float64,3}

    function Path(nbeads::Integer, nparticles::Integer; mass=1.0, spatialdimensions=3)
        R=rand(nbeads, nparticles, spatialdimensions)

        new(nbeads, nparticles, mass, spatialdimensions, R)
    end
end


function KineticAction(sA,sB)
    KE=(sA.r .- sB.r).^2
    KE=KE/4*Lambda*tau
    KE
end

function PotentialAction(sA,sB) #sliceA, sliceB
    # Primitive action
    PE=Vext(sA)+Vext(sB) + Vint(sA,sB)
    PE=PE*0.5*tau
    PE
end

