#!/bin/env python

import numpy as np

from vmc_postproc import recspace

def delta(kx,ky,params):
    return 0.5*(np.exp(1j*params['phi']*np.pi)*np.cos(kx*2*np.pi)+\
                np.exp(-1j*params['phi']*np.pi)*np.cos(ky*2*np.pi))

def mfham(kx,ky,params):
    H=np.zeros([4,4]+list(np.shape(kx)),dtype=complex)
    H[0,0,:]=-params['neel']
    H[0,1,:]=-delta(kx,ky,params).conjugate()
    H[1,0,:]=-delta(kx,ky,params)
    H[1,1,:]=params['neel']
    H[2,2,:]=params['neel']
    H[2,3,:]=-delta(kx,ky,params).conjugate()
    H[3,2,:]=-delta(kx,ky,params)
    H[3,3,:]=-params['neel']
    return H

def evals(kx,ky,params):
    v=np.zeros([4]+list(np.shape(kx)))
    v[0,:]=-np.sqrt(params['neel']**2+abs(delta(kx,ky,params))**2)
    v[1,:]=-np.sqrt(params['neel']**2+abs(delta(kx,ky,params))**2)
    v[2,:]=np.sqrt(params['neel']**2+abs(delta(kx,ky,params))**2)
    v[3,:]=np.sqrt(params['neel']**2+abs(delta(kx,ky,params))**2)
    return v

def evecs(kx,ky,params):
    v=np.zeros([4,4]+list(np.shape(kx)),dtype=complex)
    ev=evals(kx,ky,params)
    dk=delta(kx,ky,params)
    pdk=np.exp(1j*np.angle(dk))
    # spin up band -1
    v[0,0,:]=np.sqrt(0.5*(1.0-params['neel']/ev[0,:]))
    v[1,0,:]=pdk*np.sqrt(0.5*(1.0+params['neel']/ev[0,:]))
    # spin down band -1
    v[2,1,:]=np.sqrt(0.5*(1.0+params['neel']/ev[1,:]))
    v[3,1,:]=pdk*np.sqrt(0.5*(1.0-params['neel']/ev[1,:]))
    # spin up band +1
    v[0,2,:]=-pdk.conjugate()*np.sqrt(0.5*(1.0-params['neel']/ev[2,:]))
    v[1,2,:]=np.sqrt(0.5*(1.0+params['neel']/ev[2,:]))
    # spin down band +1
    v[2,3,:]=-pdk.conjugate()*np.sqrt(0.5*(1.0+params['neel']/ev[3,:]))
    v[3,3,:]=np.sqrt(0.5*(1.0-params['neel']/ev[3,:]))
    return v

def spinops(params):
    kx,ky=recspace.mbz(params)
    kqx=kx-params['qx']/params['L']
    kqy=ky-params['qy']/params['L']
    evk=evecs(kx,ky,params)
    evq=evecs(kqx,kqy,params)
    if (params['qx'],params['qy']) in [(0,0),(params['L']/2,params['L']/2)]\
            and params['channel']=='long':
        nexc=params['L']**2/2
        so=[np.zeros(nexc,dtype=complex),np.zeros(nexc,dtype=complex),np.zeros(nexc+1,dtype=complex)]
    else:
        nexc=params['L']**2
        so=[np.zeros(nexc,dtype=complex),np.zeros(nexc,dtype=complex),np.zeros(nexc,dtype=complex)]
    #Sqx=Sqpm=Sqmp=Sqy
    so[0]=np.conj(evk[0,2,:])*evq[2,1,:]+np.conj(evk[1,2,:])*evq[3,1,:]
    #Sqy=Sqx
    so[1]=so[0]
    #Sqz
    so[2][::2]=0.5*(np.conj(evk[0,2,:])*evq[0,0,:]+np.conj(evk[1,2,:])*evq[1,0,:])
    so[2][1::2]=-0.5*(np.conj(evk[2,3,:])*evq[2,1,:]+np.conj(evk[3,3,:])*evq[3,1,:])
    return so

def refstate(params):
    state=np.zeros(2*params['L']**2)
    state[::2]=1
    return state
