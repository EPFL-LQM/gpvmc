#!/bin/env python

import proc,load
import numpy as np
import re
import six

def get_eig_sys(filenames,nsamp,wav=None,tol=1e-12):
    if type(filenames)!=list:
        filenames=[filenames]
    HO=load.get_quantity(filenames,nsamp)
    H=HO[:,:HO.shape[1]/2,:]
    O=HO[:,HO.shape[1]/2:,:]
    if not wav:
        print filenames
        mo=re.search(r'(.*/[0-9]+)-ProjHeis\.h5',filenames[0])
        wav=mo.groups()[0]+'-WaveFunction.h5'
    fs=np.diag(proc.fermisigns(load.get_wav(wav)))
    H=np.einsum('ij,kjl,lm->kim',fs,H,fs)
    O=np.einsum('ij,kjl,lm->kim',fs,O,fs)
    evals=np.zeros(H.shape[:2])
    evecs=np.zeros(H.shape,dtype=complex)
    for s in range(nsamp):
        evals[s,:],evecs[s,:,:]=proc.geneigh(H[s,:,:],O[s,:,:],tol)
    return H,O,evals,evecs
