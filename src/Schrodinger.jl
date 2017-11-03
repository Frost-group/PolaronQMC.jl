# Simple discretised Schr. solver; with approx KE

# Mainly liberated from myself: https://github.com/jarvist/Julia-SoftModeTISH-DeformationPotential

function sampledV(V, N=50)
    discreteV=[ V(i/(N-1) - 1/2) for i=0:(N-1) ]
    return discreteV
end

function TISE(V,N=50,dx=1E2/(N-1))
    # PE terms on the trace
    diagonal = sampledV(V,N) .+ 2.0/dx^2
        #[(2.0/dx^2 + V(r))::Float64 for r in -range:2*range/N:range]
    # KE terms on the tridiagonals
    updiagonal = [(-1/dx^2)::Float64 for r in 1:N-1]
    
    H =diagm(diagonal,0) + diagm(updiagonal,1) + diagm(updiagonal,-1)
    
    # And so we solve; explicit diagonalisation
    evals,evecs=eig(H)
   
    return evals,evecs
end


