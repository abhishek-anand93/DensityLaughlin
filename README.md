# DensityLaughlin
This fortran code computes the density of particles in fractional quantum hall state over a circular system, defined by Laughlin state of given filling factor. Filling factor 1 gives us integer qH (quantum Hall) state and its density, where as filling factor 1/3 gives us the density of fqH state. 

'additionalpiece' and 'ap1' are just excitations in the state. ap1 is when a quasi-hole is placed at some position.

It generates a matrix where row and column have radial and angular profile of density respectively.
A python code is attached to read the file and generate a 2D heat plot.

Laughlin.f90 will take 'N' (number of particles), 'm' in the filling fraction (1/m) of the state (eg. m=3 for a ff 1/3),
'trial' to keep track of number of trials and also because Fortran gives error when it tries to generate files with same caption.



