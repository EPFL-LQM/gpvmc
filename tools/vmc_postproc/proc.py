#!/bin/env python

import copy
import numpy as np
import numpy.linalg as lg
from scipy.linalg import eigh
import warnings
from vmc_postproc import stagflux,sfpnxphz

def geneigh(A,B,tol=1e-12):
    """
    Solves the generalized eigenvalue problem also in the case where A and B share a common
    null-space. The eigenvalues corresponding to the null-space are given a Nan value.
    The null-space is defined with the tolereance tol.
    """
    # first check if there is a null-space issue
    if lg.matrix_rank(B,tol)==np.shape(B)[0]:
        return eigh(A,B)
    # first diagonalize the overlap matrix B
    Be,Bv=eigh(B)
    # rewrite the A matrix in the B-matrix eigenspace
    At=np.dot(np.conj(Bv.T),np.dot(A,Bv))
    Bt=np.diag(Be)
    # detect shared null-space. that is given by the first n null eigenvalues of B
    try:
        idx=next(i for i,v in enumerate(Be) if v>tol)
    except StopIteration:
        raise(RuntimeError('geneigh: Rank of B < B.shape[0] but null-space could not be found!'))
    # check that the B matrix null-space is shared by A.
    m=np.amax(abs(At[0:idx,:].flatten()))
    if m>tol:
        warnings.warn('Maximum non-diagonal element in A written in B null-space is bigger than the tolerance \''+str(tol)+'\'.',UserWarning)
    # diagonalize the non-null-space part of the problem
    Et,Vt=eigh(At[idx:,idx:],Bt[idx:,idx:])
    # define Ut, the change of basis in the non-truncated space
    Ut=np.zeros(np.shape(A),A.dtype)
    Ut[0:idx,0:idx]=np.eye(idx)
    Ut[idx:,idx:]=Vt
    U=np.dot(Bv,Ut)
    E=np.concatenate((float('NaN')*np.ones(idx),Et))
    return E,U

def fermisigns(states,params):
    qp_mod=None
    if params['sfpnxphz_wav']:
        qp_mod=sfpnxphz
    elif params['stagflux_wav']:
        qp_mod=stagflux
    else:
        raise RuntimeError('not implemented')
    refstate=copy.deepcopy(qp_mod.refstate(params))
    fs=np.ones(np.shape(states)[0])
    for s in range(np.shape(states)[0]):
        hole=[]
        part=[]
        ref=copy.deepcopy(refstate)
        for i in range(np.shape(states)[1]):
            if ref[i]==0 and states[s,i]!=0:
                part.append(i)
            if ref[i]!=0 and states[s,i]==0:
                hole.append(i)
            if part and hole:
                h=hole.pop()
                p=part.pop()
                fs[s]*=(-1)**(np.sum(ref[min(h,p)+1:max(h,p)]))
                ref[h]=0
                ref[p]=1
        if hole or part:
            raise(RuntimeError('fermisigns: particle number not conserved?'))
    return fs
