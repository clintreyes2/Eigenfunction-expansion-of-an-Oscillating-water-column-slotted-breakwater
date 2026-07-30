# Eigenfunction-expansion-of-an-Oscillating-water-column-slotted-breakwater
This code implements the theory described in the paper: "An analytical study of a dual-function oscillating water column with nonlinear power take-off and a slotted supporting wall" 10.1016/j.oceaneng.2026.127222. 

A non-linear PTO is performed with an iteration procedure performed until convergence of the dissipation parameter \beta_{PTO}. 

Solves the linear dispersion relation to obtain the real and imaginary wave numbers in disp_rel_ee2.m. 
Note that initial guess for the solution must be close to the right hand asymptote (closer to the positive x) and not the left hand. 

