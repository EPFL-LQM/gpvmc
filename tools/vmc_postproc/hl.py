#!/bin/env python

from vmc_postproc import proc,load,stagflux,sfpnxphz
import numpy as np
import re
import six

def get_stat_spin_struct(filenames,nsamp):
    """
    Gets the static structure factor flatened.
    The q-vectors are given by:
    q=arange(L)
    qx,qy=meshgrid(q,q)
    qx=qx.flatten()
    qy=qy.flatten()
    """
    if type(filenames)!=list:
        filenames=[filenames]
    Sq=load.get_quantity(filenames,nsamp)
    Sqxx=0.25*(Sq[:,1,:]+Sq[:,2,:]+Sq[:,3,:]+Sq[:,4,:])
    Sqyy=0.25*(Sq[:,1,:]+Sq[:,2,:]-Sq[:,3,:]-Sq[:,4,:])
    Sqzz=Sq[:,0,:]
    return Sqxx,Sqyy,Sqzz

def get_eig_sys(filenames,nsamp,wav=None,statstruct=None,tol=1e-12):
    if type(filenames)!=list:
        filenames=[filenames]
    HO=load.get_quantity(filenames,nsamp)
    H=HO[:,:HO.shape[1]/2,:]
    O=HO[:,HO.shape[1]/2:,:]
    if not wav:
        mo=re.search(r'(.*/[0-9]+)-ProjHeis\.h5',filenames[0])
        wav=mo.groups()[0]+'-WaveFunction.h5'
    wav_st=load.get_wav(wav)
    fs=np.diag(proc.fermisigns(wav_st,load.get_attr(filenames[0])))
    renorm=1.0;
    H=np.einsum('ij,kjl,lm->kim',fs,H,fs)
    O=np.einsum('ij,kjl,lm->kim',fs,O,fs)
    evals=np.zeros(H.shape[:2])
    evecs=np.zeros(H.shape,dtype=complex)
    for s in range(nsamp):
        evals[s,:],evecs[s,:,:]=proc.geneigh(H[s,:,:],O[s,:,:],tol)
    if statstruct:
        params=load.get_attr(filenames[0])
        qs=[]
        for qx in range(int(params['L'])):
            for qy in range(int(params['L'])):
                qs.append((qx,qy))
        idx=qs.index((int(params['qx']),int(params['qy'])))
        if params['stagflux_wav']:
            pkq=stagflux.spinops(params)
            if params['channel']=='trans':
                sqpm=np.einsum('i,sij,j->s',np.conj(pkq[1]),O,pkq[1])
                for s in range(nsamp):
                    O[s,:,:]*=statstruct[1][0,idx]*2/sqpm
            elif params['channel']=='long':
                sqzz=np.einsum('i,sij,j->s',np.conj(pkq[2]),O,pkq[2])
                for s in range(nsamp):
                    O[s,:,:]*=statstruct[2][0,idx]*2/sqzz
            else:
                raise RuntimeError('Channel not undertood')
        elif params['sfpnxphz_wav']:
            pkq=sfpnxphz.spinops(params)
            sqzz=np.einsum('i,sij,j->s',np.conj(pkq[2]),O,pkq[2])
            for s in range(nsamp):
                O[s,:,:]*=statstruct[2][0,idx]*2/sqzz
        else:
            raise RuntimeError('Not implemented')
    return H,O,evals,evecs

def get_sq_ampl(filenames,nsamp,wav=None,tol=1e-12,H=None,O=None,evals=None,evecs=None):
    if type(filenames)!=list:
        filenames=[filenames]
    if H==None:
        H,O,evals,evecs=get_eig_sys(filenames,nsamp,wav,tol)
    qp_mod=None
    params=load.get_attr(filenames[0])
    if params['sfpnxphz_wav']:
        qp_mod=sfpnxphz
    elif params['stagflux_wav']:
        qp_mod=stagflux
    else:
        raise RuntimeError('sfpnzphx_wav not implemented yet')
    pkq=qp_mod.spinops(params)
    if params['stagflux_wav']:
        if params['channel']=='trans':
            return abs(np.einsum('i,sij,sjk->sk',np.conj(pkq[1]),O,evecs))**2
        elif params['channel']=='long':
            return abs(np.einsum('i,sij,sjk->sk',np.conj(pkq[2]),O,evecs))**2
        else:
            raise RuntimeError('channel not understood.')
    elif params['sfpnxphz_wav']:
        out=[None]*3
        for a in range(3):
            out[a]=abs(np.einsum('i,sij,sjk->sk',np.conj(pkq[a]),O,evecs))**2
        return out
