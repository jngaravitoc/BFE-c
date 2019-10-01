BFE-c
==========

The code computes coefficients of the [Hernquist Basis Field Expansion](https://ui.adsabs.harvard.edu/abs/1992ApJ...386..375H/abstract)  from snapshots of N-body simulations. 
At the moment, the code only reads outputs from Gadget. 

BFE-c is wirtten in C and it can run in parallel for an efficient computation of the coefficients.
In addition, the code also computes the variance matrix of the coefficients (following Weinberg+98) 
which will be useful for post noise-reduction post-processing. 

Noise-reduction:
---------------
BFE-c randomly samples the particle 
distribution and compute the coefficients in each of the random realizations. The final output is the 
mean of the coefficients in each of the realizations. This procedure minimizes the noise
from the discrete particle distribution. 


Install:
--------


Dependencies:
-------------

- gsl
- openmp 

How to run:
-----------



