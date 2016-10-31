# HALCYON - PIGS might fly

Writing a Path Integral Ground State monte-carlo code in the Julia programming language. A work in progress!

## References

Textbook: "Interacting Electrons: Theory and Computational Approaches" 
Richard M. Martin, Lucia Reining, and David M. Ceperley. 
Cambridge University Press. (2016)
Section 25.4 (Ground-state path integrals (GSPI))
http://dx.doi.org/10.1017/CBO9781139050807 

PIGS papers:

[1] "Reptation Quantum Monte Carlo: A Method for Unbiased Ground-State Averages and
Imaginary-Time Correlations"
Stefano Baroni and Saverio Moroni
Phys. Rev. Lett. 82, 4745 (1999); 
http://dx.doi.org/10.1103/PhysRevLett.82.4745

[2] "Path integrals in the theory of condensed helium" 
D. M. Ceperley
Rev. Mod. Phys. 67, 279 (1995); 
http://dx.doi.org/10.1103/RevModPhys.67.279

[3] A path integral ground state method
A. Sarsa, K. E. Schmidt and W. R. Magro
J. Chem. Phys. 113, 1366 (2000); 
http://dx.doi.org/10.1063/1.481926

Fermion sign problem papers:

[1] Sign problem in the numerical simulation of many-electron systems
E. Y. Loh, Jr., J. E. Gubernatis, R. T. Scalettar, S. R. White, D. J. Scalapino, and R. L. Sugar
Phys. Rev. B 41, 9301 – Published 1 May 1990
https://doi.org/10.1103/PhysRevB.41.9301

[2] Monte Carlo simulations with indefinite and complex-valued measures
T. D. Kieu and C. J. Griffin
Phys. Rev. E 49, 3855 – Published 1 May 1994
https://doi.org/10.1103/PhysRevE.49.3855

## Lectures and other notes (including finite T PIMC)

2012 - quantum Monte Carlo: Theory and Fundies - Uni of illinois
http://www.mcc.uiuc.edu/summerschool/2012/program.html

BryanClark: Tutorial: Writing a Path Integral Code in Python
http://web.engr.illinois.edu/~bkclark/PIMCTutorial/tutorial.pdf

Adrian Del Maestro - 14 pages of notes and a beautifully clean 222 line PIMC Python code for a 1D SHO 
https://github.com/agdelma/pimc-notes

## Other codes

See also, @amaciarey's F90 PIGS code, for dipolar bosons in 2D: 
https://github.com/amaciarey/PathIntegralGroundState

PIMC++: finite T PIMC code, C++: 
https://github.com/bkclark/pimcpp

ALSO; same author (?); Simple PIMC
https://github.com/etano/simpimc
