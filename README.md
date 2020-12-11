# DensityLaughlin
This fortran code computes the density of particles in fractional quantum hall state in a disk geometry, defined by Laughlin state of given filling factor =1/m. For instance, filling factor 1 generates density of integer qH (quantum Hall)

'additionalpiece' and 'ap1' are just excitations in the state. ap1 is when a quasi-hole is placed at some position.

It generates a matrix where row and column have radial and angular profile of density respectively.
A python code is attached to read the file and generate a 2D heat plot.

Laughlin.f90 will take 'N' (number of particles), 'm' in the filling fraction (1/m) of the state (eg. m=3 for a ff 1/3),
'trial' to keep track of number of trials (averaging over trials with different seeds will smoothen the density faster!).

angular_plot.py will take N, m and 'trial' as input and generate the plot.
