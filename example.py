#/bin/env python

helpstr="""
# This is a script showing an example to obtain
# the ground state energy, staggered magnetizatio
# and the static spin structure factor of the |SF+N>
# wavefunction.
# In a second step the q=(pi,0) component of the
# dynamical spin structure factor S(q,w) is calculated.
# This script requires numpy, matplotlib, six and h5py installed.
#
# The postprocessing tools coming with the variational
# monte carlo source code must be installed in a location
# known to python.
#
# The \"vmc\" executable must be installed in a location in the PATH
# environment variable. It must have been built with the USEMPI=NO
# option.
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import re
import six

from vmc_postproc import hl,load,proc

print(helpstr)

# ground state calculation

command="vmc --neel=0.055 --phi=0.085 --L=8 --samples=1 --samples_saves=1 --samples_saves_stat=10000 --therm=200 --channel=groundstate --meas_projheis --meas_stagmagn --meas_statspinstruc".split(' ')

pgs=subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
stdout,stderr=pgs.communicate()
if six.PY3:
    stdout=str(stdout,encoding='utf-8')

# obtain groundstate calculation output files prefix

gsprefix=re.search(r'output prefix=([0-9]+)',stdout).groups()[0]

gsen,_,_,_=hl.get_eig_sys('{}-ProjHeis.h5'.format(gsprefix),1)
stagmagn=load.get_quantity('{}-StagMagnZ.h5'.format(gsprefix),1)
sqaa,sraa=hl.get_stat_spin_struct('{}-StatSpinStruct.h5'.format(gsprefix),1)

print("""Ground state energy: {gsen},
Ground state staggered magnetization: {stz},
Plotting <Sx(0)Sx(r)> correlation function in sqxx.png""".format(gsen=np.real(gsen[0,0,0]),stz=np.real(stagmagn[0,0,0])))

plt.figure()
plt.pcolormesh(np.arange(-4,4)-0.5,np.arange(-4,4)-0.5,np.real(sraa[0][0,:,:]))
plt.colorbar()
plt.savefig('srxx.png')

# excited state calculation for momentum q=(pi,0)

command="vmc --neel=0.055 --phi=0.085 --L=8 --samples=1 --samples_saves=1 --samples_saves_stat=10000 --therm=200 --channel=trans --meas_projheis --qx=4 --qy=0".split()
pgs=subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
stdout,stderr=pgs.communicate()
if six.PY3:
    stdout=str(stdout,encoding='utf-8')

# obtain excited state calculation output file prefix

excprefix=re.search(r'output prefix=([0-9]+)',stdout).groups()[0]

H,O,evals,evecs=hl.get_eig_sys('{}-ProjHeis.h5'.format(excprefix),1,wav='{}-WaveFunction.h5'.format(excprefix),statstruct=sqaa)

sq=hl.get_sq_ampl('{}-ProjHeis.h5'.format(excprefix),1,wav='{}-WaveFunction.h5'.format(excprefix),O=O,evecs=evecs)

# convolute the excitations amplitudes with gaussians for better visualization

w=np.arange(0,6,0.001) # in units of J

sqw=proc.gaussians(w,np.real(evals[0,:])-64*np.real(gsen[0,0,0]),sq[0,:],0.1)

print('Plotting transverse S(q=(pi,0),w) in sqwx.png')

plt.figure()
plt.plot(w,sqw,label=r'$q=(\pi,0)$')
plt.setp(plt.gca(),xlabel=r'$\omega$ [units of J]',ylabel=r'$S(q,\omega)$')
plt.savefig('sqwx.png')
