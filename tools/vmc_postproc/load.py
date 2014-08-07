#!/bin/env python
from __future__ import unicode_literals,print_function
import six
import subprocess
import sys
import re
import h5py
import copy
import numpy as np
import numpy.linalg as lg
from scipy.linalg import eigh
import warnings

def get_quantity(filepath,Nsamp):
    if type(filepath)!=list:
        fielpath=[filepath]
    dpath,args=get_stat(filepath,Nsamp)
    hfile=h5py.File(dpath[args[0][0]][0],'r')
    dshape=list(hfile[dpath[args[0][0]][1]].shape)
    hfile.close()
    dshape[1]/=2 #because complex numbers a contiguous along a line
    Q=np.zeros([Nsamp]+list(dshape),dtype=complex)
    for si,bi in enumerate(args):
        for d in bi:
            hfile=h5py.File(dpath[d][0],'r')
            dat=hfile[dpath[d][1]]
            Q[si,:]+=dat[:,::2]+1j*dat[:,1::2]
            hfile.close()
        Q[si,:]/=len(bi)
    return Q

def get_stat(filename,Nsamp=1):
    hfile=[]
    if type(filename)==list:
        for i,f in enumerate(filename):
            hfile.append(h5py.File(f,'r'))
    elif type(filename)==str:
        hfile.append(h5py.File(filename,'r'))
        filename=[filename]
    stats=[]
    datapath=[]
    for ih,h in enumerate(hfile):
        for r in h:
            for d in h[r]:
                try:
                    stats.append(int(h[r][d].attrs['statistics'][0]))
                except KeyError as err:
                    stats.append(1)
                datapath.append((filename[ih],"/{0}/{1}".format(r,d)))
    bunches,args=bunch(stats,Nsamp,indices=True)
    addstat=np.array([sum(bunches[i]) for i in range(Nsamp)])
    print("Average statistics of {0} (min: {1}, max: {2})"\
            .format(np.mean(addstat),np.amin(addstat),np.amax(addstat)))
    for f in hfile:
        f.close()
    return datapath,args

def bunch(instat,Nsamp,indices=False):
    """
    Bunches the integers from instat into Nsamp bunches
    while trying to even out the sum of the bunches integers.
    This is not garanteed to be the best solution. It is just
    a first guess, may probably be optimized further...
    Returns the bunched integers.
    if indices==True, also returns the bunched indices corresponding
    to the original instat list.
    """
    stats=list(instat)
    args=argsort(stats)
    stats.sort()
    bunched=[[] for n in range(Nsamp)]
    argbunched=copy.deepcopy(bunched)
    n=0
    while len(stats)>2*Nsamp-1:
        bunched[n]+=[stats[0],stats[-1]]
        argbunched[n]+=[args[0],args[-1]]
        stats.pop()
        stats.pop(0)
        args.pop()
        args.pop(0)
        n+=1
        if n==Nsamp:
            n=0
    bunchedstat=[sum(b) for b in bunched]
    while len(stats)>0:
        bunched[np.argmin(bunchedstat)]+=[stats[-1]]
        argbunched[np.argmin(bunchedstat)]+=[args[-1]]
        bunchedstat[np.argmin(bunchedstat)]+=stats[-1]
        stats.pop()
        args.pop()
    if not indices:
        return bunched
    else:
        return bunched,argbunched

def argsort(seq):
    return sorted(range(len(seq)),key=seq.__getitem__)

def get_attr(filename):
    hfile=h5py.File(filename,'r')
    attrs=dict(hfile.attrs)
    hfile.close()
    if six.PY3:
        for k,v in attrs.items():
            val=v
            if type(v) is np.bytes_:
                val=str(val,encoding='ascii')
            attrs[k]=val
    return attrs

def get_wav(filename):
    attr=get_attr(filename)
    filetype=''
    try:
        filetype=attr['type']
    except KeyError:
        mo=re.match('.*/?[0-9]+-(.*)\.h5',filename)
        filetype=mo.groups()[0]
    if filetype!='WaveFunction':
        raise(RuntimeError('get_wav: Bad file type \'{}\'.'.format(filetype)))
    L=attr['L']
    hfile=h5py.File(filename,'r')
    states_flav=[]
    for k in hfile.keys():
        states_flav+=[hfile[k]]
    states=np.column_stack(states_flav)
    hfile.close()
    return states
