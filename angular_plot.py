#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
from os import sys
import matplotlib as mpl
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
#=========================================

N=int(sys.argv[1]);
filling=int(sys.argv[2]);
trial=int(sys.argv[3]);

#=====================================================================

filename = 'density_N%d_filling%d_trial%d'	% (N,filling,trial);
filename2 = 'init_%d_%d_%d'                     % (N,filling,trial);

if exists(filename2):
    a=np.loadtxt(filename2)
    NUMBINS=int(np.round(a[0]));
    print(NUMBINS)
    dR2=a[1]
    dtheta=a[2]
    rmax = a[3];

#===========================================================================
'''This part reads the 2D array of NUMBIN*NUMBIN, row by row.
   Thus, the returned list y consists of a linear array of dimension
   [NUMBIN,0] with each element is itself a list containing each row.'''

def read_file(filename):
    with open(filename,'r') as data:
        x = []
        y = []
        for line in data:
            p = line.split()
            for i in range(NUMBINS):
                x.append(float(p[i]))
            y.append(np.array(x))
            x=[]
    return y

density = np.array(read_file(filename))
#==========================================================================

density = density*NUMBINS/(np.pi*dR2)
localfillingfraction=2*np.pi*density;
r = np.sqrt(np.arange(0,rmax*rmax,dR2) )[:NUMBINS]
theta = np.linspace(0,2*np.pi,NUMBINS)

mpl.rcParams['legend.fontsize'] = 10

#===============================================================
'''the matrix needed to be trasposed, to get the radial and angular
   data to be read correctly'''

sdensity = np.transpose(localfillingfraction)

fig = plt.figure()
ax = Axes3D(fig)
rad, th = np.meshgrid(r,theta)
plt.subplot(projection="polar")
plt.pcolormesh(th,rad,sdensity)
plt.thetagrids([theta * 15 for theta in range(360//15)])
plt.plot(theta,rad,color='k',ls='none')
plt.grid()
plt.savefig("N%d_FILLING%d_trial%d.pdf" % (N,filling,trial));
plt.show()
