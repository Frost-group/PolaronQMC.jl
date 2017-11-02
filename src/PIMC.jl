# PIMC.jl
# Initially location for all the path integral bits and bobs

# Very generic path
struct Path
    timeslice
    particle
    spatialdimension
end

const T=300
const tau=0.5

timeslices=round(1/(tau*T))
particles=2
spatialdimensions=1

path=zeros(timeslices,particles,spatialdimensions)
show(path)

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
