#!/bin/env python

import numpy as np
from vmc_postproc import recspace

def delta(kx,ky,params):
    return 0.5*(np.exp(1j*params['phi'])*np.cos(kx*2*np.pi)+\
                np.exp(-1j*params['phi'])*np.cos(ky*2*np.pi))

def mfham(kx,ky,params):
    H=np.zeros([4,4]+list(np.shape(kx)),dtype=complex)
    dk=delta(kx,ky,params)
    H[0,0,:]=params['field']
    H[0,1,:]=-dk.conjugate()
    H[0,2,:]=params['neel']
    H[1,0,:]=-dk
    H[1,1,:]=params['field']
    H[1,3,:]=-params['neel']
    H[2,0,:]=params['neel']
    H[2,2,:]=-params['field']
    H[2,3,:]=-dk.conjugate()
    H[3,1,:]=-params['neel']
    H[3,2,:]=-dk
    H[3,3,:]=-params['field']
    Tp=trans(kx,ky,params)
    return np.einsum('ij...,jk...,lk...->il...',Tp,H,Tp)

def evals(kx,ky,params):
    v=np.zeros([4]+list(np.shape(kx)))
    v[0,:]=-np.sqrt(params['neel']**2+(params['field']-abs(delta(kx,ky,params)))**2)
    v[1,:]=-np.sqrt(params['neel']**2+(params['field']+abs(delta(kx,ky,params)))**2)
    v[2,:]=np.sqrt(params['neel']**2+(params['field']-abs(delta(kx,ky,params)))**2)
    v[3,:]=np.sqrt(params['neel']**2+(params['field']+abs(delta(kx,ky,params)))**2)
    return v

def evecs(kx,ky,params):
    v=np.zeros([4,4]+list(np.shape(kx)),dtype=complex)
    ev=evals(kx,ky,params)
    dk=delta(kx,ky,params)
    # band 1
    v[0,0,:]=0.5*np.sqrt(1.0-params['neel']/evals(kx,ky,params)[0,:])
    v[1,0,:]=0.5*np.sign(abs(dk)-params['field'])*np.exp(1j*np.angle(dk))*np.sqrt(1.0+params['neel']/evals(kx,ky,params)[0,:])
    v[2,0,:]=0.5*np.sign(abs(dk)-params['field'])*np.sqrt(1.0+params['neel']/evals(kx,ky,params)[0,:])
    v[3,0,:]=0.5*np.exp(1j*np.angle(dk))*np.sqrt(1.0-params['neel']/evals(kx,ky,params)[0,:])
    # band 2
    v[0,1,:]=0.5*np.sqrt(1.0-params['neel']/evals(kx,ky,params)[1,:])
    v[1,1,:]=0.5*np.sign(abs(dk)+params['field'])*np.exp(1j*np.angle(dk))*np.sqrt(1.0+params['neel']/evals(kx,ky,params)[1,:])
    v[2,1,:]=-0.5*np.sign(abs(dk)+params['field'])*np.sqrt(1.0+params['neel']/evals(kx,ky,params)[1,:])
    v[3,1,:]=-0.5*np.exp(1j*np.angle(dk))*np.sqrt(1.0-params['neel']/evals(kx,ky,params)[1,:])
    # band 3
    v[0,2,:]=0.5*np.sqrt(1.0-params['neel']/evals(kx,ky,params)[2,:])
    v[1,2,:]=-0.5*np.sign(abs(dk)-params['field'])*np.exp(1j*np.angle(dk))*np.sqrt(1.0+params['neel']/evals(kx,ky,params)[2,:])
    v[2,2,:]=-0.5*np.sign(abs(dk)-params['field'])*np.sqrt(1.0+params['neel']/evals(kx,ky,params)[2,:])
    v[3,2,:]=0.5*np.exp(1j*np.angle(dk))*np.sqrt(1.0-params['neel']/evals(kx,ky,params)[2,:])
    # band 4
    v[0,3,:]=0.5*np.sqrt(1.0-params['neel']/evals(kx,ky,params)[3,:])
    v[1,3,:]=-0.5*np.sign(params['field']+abs(dk))*np.exp(1j*np.angle(dk))*np.sqrt(1.0+params['neel']/evals(kx,ky,params)[3,:])
    v[2,3,:]=0.5*np.sign(params['field']+abs(dk))*np.sqrt(1.0+params['neel']/evals(kx,ky,params)[3,:])
    v[3,3,:]=-0.5*np.exp(1j*np.angle(dk))*np.sqrt(1.0-params['neel']/evals(kx,ky,params)[3,:])
    return v
    #return np.einsum('ij...->ji...',v).conjugate()

def trans(kx,ky,params):
    Tp=np.zeros([4,4]+list(np.shape(kx)))
    Tp[0,0,:]=1
    Tp[0,2,:]=-1
    Tp[1,1,:]=1
    Tp[1,3,:]=-1
    Tp[2,0,:]=1
    Tp[2,2,:]=1
    Tp[3,1,:]=1
    Tp[3,3,:]=1
    return Tp/np.sqrt(2)

def spinops(params):
    kx,ky=recspace.mbz(params)
    kqx,kqy=recspace.mbzmod(kx-params['qx']/params['L'],ky-params['qy']/params['L'])
    skq=2*np.array(recspace.inmbz(kx-params['qx']/params['L'],ky-params['qy']/params['L']),dtype=int)-1
    evk=evecs(kx,ky,params)
    evq=evecs(kqx,kqy,params)
    nexc=2*params['L']**2
    s=0
    if (params['qx'],params['qy']) in [(0,0),(params['L']/2,params['L']/2)]:
        nexc+=1
        s=1
    so=[np.zeros((nexc),dtype=complex)]*3
    #Sqx
    so[0][s::4]=0.5*(-evk[0,2,:].conjugate()*evq[0,0,:]+evk[2,2,:].conjugate()*evq[2,0,:]+\
                 skq*(-evk[1,2,:].conjugate()*evq[1,0,:]+evk[3,2,:].conjugate()*evq[3,0,:]))
    so[0][s+1::4]=0.5*(-evk[0,3,:].conjugate()*evq[0,0,:]+evk[2,3,:].conjugate()*evq[2,0,:]+\
                 skq*(-evk[1,3,:].conjugate()*evq[1,0,:]+evk[3,3,:].conjugate()*evq[3,0,:]))
    so[0][s+2::4]=0.5*(-evk[0,2,:].conjugate()*evq[0,1,:]+evk[2,2,:].conjugate()*evq[2,1,:]+\
                 skq*(-evk[1,2,:].conjugate()*evq[1,1,:]+evk[3,2,:].conjugate()*evq[3,1,:]))
    so[0][s+3::4]=0.5*(-evk[0,3,:].conjugate()*evq[0,1,:]+evk[2,3,:].conjugate()*evq[2,1,:]+\
                 skq*(-evk[1,3,:].conjugate()*evq[1,1,:]+evk[3,3,:].conjugate()*evq[3,1,:]))
    #Sqy
    so[1][s::4]=0.5j*(evk[0,2,:].conjugate()*evq[2,0,:]-evk[2,2,:].conjugate()*evq[0,0,:]+\
                skq*(evk[1,2,:].conjugate()*evq[3,0,:]-evk[3,2,:].conjugate()*evq[1,0,:]))
    so[1][s+1::4]=0.5j*(evk[0,3,:].conjugate()*evq[2,0,:]-evk[2,3,:].conjugate()*evq[0,0,:]+\
                  skq*(evk[1,3,:].conjugate()*evq[3,0,:]-evk[3,3,:].conjugate()*evq[1,0,:]))
    so[1][s+2::4]=0.5j*(evk[0,2,:].conjugate()*evq[2,1,:]-evk[2,2,:].conjugate()*evq[0,1,:]+\
                  skq*(evk[1,2,:].conjugate()*evq[3,1,:]-evk[3,2,:].conjugate()*evq[1,1,:]))
    so[1][s+3::4]=0.5j*(evk[0,3,:].conjugate()*evq[2,1,:]-evk[2,3,:].conjugate()*evq[0,1,:]+\
                  skq*(evk[1,3,:].conjugate()*evq[3,1,:]-evk[3,3,:].conjugate()*evq[1,1,:]))
    #Sqz
    so[2][s::4]=0.5*(evk[0,2,:].conjugate()*evq[2,0,:]+evk[2,2,:].conjugate()*evq[0,0,:]+\
               skq*(evk[1,2,:].conjugate()*evq[3,0,:]+evk[3,2,:].conjugate()*evq[1,0,:]))
    so[2][s+1::4]=0.5*(evk[0,3,:].conjugate()*evq[2,0,:]+evk[2,3,:].conjugate()*evq[0,0,:]+\
                 skq*(evk[1,3,:].conjugate()*evq[3,0,:]+evk[3,3,:].conjugate()*evq[1,0,:]))
    so[2][s+2::4]=0.5*(evk[0,2,:].conjugate()*evq[2,1,:]+evk[2,2,:].conjugate()*evq[0,1,:]+\
                 skq*(evk[1,2,:].conjugate()*evq[3,1,:]+evk[3,2,:].conjugate()*evq[1,1,:]))
    so[2][s+3::4]=0.5*(evk[0,3,:].conjugate()*evq[2,1,:]+evk[2,3,:].conjugate()*evq[0,1,:]+\
                 skq*(evk[1,3,:].conjugate()*evq[3,1,:]+evk[3,3,:].conjugate()*evq[1,1,:]))
    if (params['qx'],params['qy']) in [(0,0),(params['L']/2,params['L']/2)]:
        so[0][0]=0.5*np.sum(-evk[0,0,:].conjugate()*evq[0,0,:]+evk[2,0,:].conjugate()*evq[2,0,:]+\
                    skq*(-evk[1,0,:].conjugate()*evq[1,0,:]+evk[3,0,:].conjugate()*evq[3,0,:]))+\
                 0.5*np.sum(-evk[0,1,:].conjugate()*evq[0,1,:]+evk[2,1,:].conjugate()*evq[2,1,:]+\
                    skq*(-evk[1,1,:].conjugate()*evq[1,1,:]+evk[3,1,:].conjugate()*evq[3,1,:]))
    return so

def refstate(params):
    state=np.zeros(2*params['L']**2)
    state[::4]=1
    state[1::4]=1
    return state
