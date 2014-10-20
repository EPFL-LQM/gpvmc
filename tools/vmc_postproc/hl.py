#!/bin/env python

from vmc_postproc import proc,load,stagflux,sfpnxphz
import numpy as np
import numpy.fft as fft
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
    params=load.get_attr(filenames[0])
    Lx=int(params['L'])
    Ly=int(params['L'])
    hLx=int(Lx/2)
    hLy=int(Ly/2)
    N=Lx*Ly
    if Sq.shape[2]==N:
        # old file format, struct stored in Fourier components
        Sqxx=np.reshape(0.25*(Sq[:,1,:]+Sq[:,2,:]+Sq[:,3,:]+Sq[:,4,:]),(Sq.shape[0],Lx,Ly))
        Sqyy=np.reshape(0.25*(Sq[:,1,:]+Sq[:,2,:]-Sq[:,3,:]-Sq[:,4,:]),(Sq.shape[0],Lx,Ly))
        Sqzz=np.reshape(Sq[:,0,:],(Sq.shape[0],Lx.Ly))
        Srxx=fft.fftshift(fft.fft2(Sqxx,axes=(1,2)),axes=(1,2))/N
        Sryy=fft.fftshift(fft.fft2(Sqyy,axes=(1,2)),axes=(1,2))/N
        Srzz=fft.fftshift(fft.fft2(Sqzz,axes=(1,2)),axes=(1,2))/N
    else :
        # new file format, struct stored as real space site pairs.
        rx,ry=np.meshgrid(np.arange(Lx,dtype=int),np.arange(Ly,dtype=int))
        rx=rx.ravel()
        ry=ry.ravel()
        rix,rjx=np.meshgrid(rx,rx)
        riy,rjy=np.meshgrid(ry,ry)
        rijx=rjx-rix
        rijy=rjy-riy
        rijx[rijx>=hLx]-=Lx
        rijx[rijx<-hLx]+=Lx
        rijy[rijy>=hLy]-=Ly
        rijy[rijy<-hLy]+=Ly
        rijx=rijx.ravel()
        rijy=rijy.ravel()
        Sr=np.zeros((Sq.shape[0],5,N),dtype=complex)
        for samp in range(Sq.shape[0]):
            for t in range(N):
                Sr[samp,:,t]=np.sum(Sq[samp,:,np.where((rijy+hLy)*Lx+rijx+hLx==t)[0]],axis=0)/N
        Srxx=np.zeros((Sq.shape[0],Lx,Ly),dtype=complex)
        Sryy=np.zeros((Sq.shape[0],Lx,Ly),dtype=complex)
        Srzz=np.zeros((Sq.shape[0],Lx,Ly),dtype=complex)
        for samp in range(Sq.shape[0]):
            Srxx[samp,:,:]=np.reshape(0.25*np.sum(Sr[samp,1:,:],axis=0),(Lx,Ly))
            Sryy[samp,:,:]=np.reshape(0.25*(np.sum(Sr[samp,1:3,:],axis=0)-np.sum(Sr[samp,3:,:],axis=0)),(Lx,Ly))
            Srzz[samp,:,:]=np.reshape(Sr[samp,0,:],(Lx,Ly))
        Sqxx=fft.ifft2(fft.fftshift(Srxx,axes=(1,2)),axes=(1,2))*np.sqrt(N)
        Sqyy=fft.ifft2(fft.fftshift(Sryy,axes=(1,2)),axes=(1,2))*np.sqrt(N)
        Sqzz=fft.ifft2(fft.fftshift(Srzz,axes=(1,2)),axes=(1,2))*np.sqrt(N)
    return (Sqxx,Sqyy,Sqzz),(Srxx,Sryy,Srzz)

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
    for k in range(nsamp):
        H[k,:,:]=np.dot(fs,np.dot(H[k,:,:],fs))
        O[k,:,:]=np.dot(fs,np.dot(O[k,:,:],fs))
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
        ssqx=statstruct[0]
        ssqz=statstruct[2]
        if len(ssqx.shape) != 2:
            ssqx=ssqx[0,...].ravel()
            ssqz=ssqz[0,...].ravel()
        if params['stagflux_wav']:
            pkq=stagflux.spinops(params)
            if params['channel']=='trans':
                sqpm=np.einsum('i,sij,j->s',np.conj(pkq[1]),O,pkq[1])
                for s in range(nsamp):
                    O[s,:,:]*=ssqx[idx]*2/sqpm
            elif params['channel']=='long':
                sqzz=np.einsum('i,sij,j->s',np.conj(pkq[2]),O,pkq[2])
                for s in range(nsamp):
                    O[s,:,:]*=ssqz[idx]*2/sqzz
            else:
                raise RuntimeError('Channel not undertood')
        elif params['sfpnxphz_wav']:
            pkq=sfpnxphz.spinops(params)
            sqzz=np.einsum('i,sij,j->s',np.conj(pkq[2]),O,pkq[2])
            for s in range(nsamp):
                O[s,:,:]*=ssqz[idx]*2/sqzz
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
            out[a]=np.zeros((nsamp,O.shape[1]))
            for k in range(nsamp):
                out[a][k,:]=abs(np.dot(np.conj(pkq[a]),np.dot(O[k,:,:],evecs[k,:,:])))**2
        return out
