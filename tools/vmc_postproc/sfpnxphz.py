#!/bin/env python

import numpy as np
import recspace

def delta(kx,ky,params):
    return 0.5*(np.exp(1j*params['phi']*np.pi)*np.cos(kx*2*np.pi)+\
                np.exp(-1j*params['phi']*np.pi)*np.cos(ky*2*np.pi))

def evals(kx,ky,params):
    v=np.zeros([4]+list(np.shape(kx)))
    v[0,:]=-np.sqrt(params['nx']**2+(params['hz']-abs(delta(kx,ky,params)))**2)
    v[1,:]=-np.sqrt(params['nx']**2+(params['hz']+abs(delta(kx,ky,params)))**2)
    v[2,:]=np.sqrt(params['nx']**2+(params['hz']-abs(delta(kx,ky,params)))**2)
    v[3,:]=np.sqrt(params['nx']**2+(params['hz']+abs(delta(kx,ky,params)))**2)
    return v

def evecs(kx,ky,params):
    v=np.zeros([4,4]+list(np.shape(kx)),dtype=complex)
    ev=evals(kx,ky,params)
    dk=delta(kx,ky,params)
    # band 1
    v[0,0,:]=0.5*np.sqrt(1.0-params['nx']/evals(kx,ky,params)[0,:])
    v[1,0,:]=0.5*np.sign(abs(dk)-params['hz'])*np.exp(1j*np.angle(dk))*np.sqrt(1.0+params['nx']/evals(kx,ky,params)[0,:])
    v[2,0,:]=0.5*np.sign(abs(dk)-params['hz'])*np.sqrt(1.0+params['nx']/evals(kx,ky,params)[0,:])
    v[3,0,:]=0.5*np.exp(1j*np.angle(dk))*np.sqrt(1.0-params['nx']/evals(kx,ky,params)[0,:])
    # band 2
    v[0,1,:]=0.5*np.sqrt(1.0-params['nx']/evals(kx,ky,params)[1,:])
    v[1,1,:]=0.5*np.sign(abs(dk)+params['hz'])*np.exp(1j*np.angle(dk))*np.sqrt(1.0+params['nx']/evals(kx,ky,params)[1,:])
    v[2,1,:]=-0.5*np.sign(abs(dk)+params['hz'])*np.sqrt(1.0+params['nx']/evals(kx,ky,params)[1,:])
    v[3,1,:]=-0.5*np.exp(1j*np.angle(dk))*np.sqrt(1.0-params['nx']/evals(kx,ky,params)[1,:])
    # band 3
    v[0,2,:]=0.5*np.sqrt(1.0-params['nx']/evals(kx,ky,params)[2,:])
    v[1,2,:]=-0.5*np.sign(abs(dk)-params['hz'])*np.exp(1j*np.angle(dk))*np.sqrt(1.0+params['nx']/evals(kx,ky,params)[2,:])
    v[2,2,:]=-0.5*np.sign(abs(dk)-params['hz'])*np.sqrt(1.0+params['nx']/evals(kx,ky,params)[2,:])
    v[3,2,:]=0.5*np.exp(1j*np.angle(dk))*np.sqrt(1.0-params['nx']/evals(kx,ky,params)[2,:])
    # band 4
    v[0,3,:]=0.5*np.sqrt(1.0-params['nx']/evals(kx,ky,params)[3,:])
    v[1,3,:]=-0.5*np.sign(params['hz']+abs(dk))*np.exp(1j*np.angle(dk))*np.sqrt(1.0+params['nx']/evals(kx,ky,params)[3,:])
    v[2,3,:]=0.5*np.sign(params['hz']+abs(dk))*np.sqrt(1.0+params['nx']/evals(kx,ky,params)[3,:])
    v[3,3,:]=-0.5*np.exp(1j*np.angle(dk))*np.sqrt(1.0-params['nx']/evals(kx,ky,params)[3,:])
    return v
    #return np.einsum('ij...->ji...',v).conjugate()

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
    so=np.zeros((nexc,3),dtype=complex)
    #Sqx
    so[s::4,0]=0.5*(-evk[0,2,:].conjugate()*evq[0,0,:]+evk[2,2,:].conjugate()*evq[2,0,:]+\
               skq*(-evk[1,2,:].conjugate()*evq[1,0,:]+evk[3,2,:].conjugate()*evq[3,0,:]))
    so[s+1::4,0]=0.5*(-evk[0,3,:].conjugate()*evq[0,0,:]+evk[2,3,:].conjugate()*evq[2,0,:]+\
                 skq*(-evk[1,3,:].conjugate()*evq[1,0,:]+evk[3,3,:].conjugate()*evq[3,0,:]))
    so[s+2::4,0]=0.5*(-evk[0,2,:].conjugate()*evq[0,1,:]+evk[2,2,:].conjugate()*evq[2,1,:]+\
                 skq*(-evk[1,2,:].conjugate()*evq[1,1,:]+evk[3,2,:].conjugate()*evq[3,1,:]))
    so[s+3::4,0]=0.5*(-evk[0,3,:].conjugate()*evq[0,1,:]+evk[2,3,:].conjugate()*evq[2,1,:]+\
                 skq*(-evk[1,3,:].conjugate()*evq[1,1,:]+evk[3,3,:].conjugate()*evq[3,1,:]))
    #Sqy
    so[s::4,1]=0.5j*(evk[0,2,:].conjugate()*evq[2,0,:]-evk[2,2,:].conjugate()*evq[0,0,:]+\
                skq*(evk[1,2,:].conjugate()*evq[3,0,:]-evk[3,2,:].conjugate()*evq[1,0,:]))
    so[s+1::4,1]=0.5j*(evk[0,3,:].conjugate()*evq[2,0,:]-evk[2,3,:].conjugate()*evq[0,0,:]+\
                  skq*(evk[1,3,:].conjugate()*evq[3,0,:]-evk[3,3,:].conjugate()*evq[1,0,:]))
    so[s+2::4,1]=0.5j*(evk[0,2,:].conjugate()*evq[2,1,:]-evk[2,2,:].conjugate()*evq[0,1,:]+\
                  skq*(evk[1,2,:].conjugate()*evq[3,1,:]-evk[3,2,:].conjugate()*evq[1,1,:]))
    so[s+3::4,1]=0.5j*(evk[0,3,:].conjugate()*evq[2,1,:]-evk[2,3,:].conjugate()*evq[0,1,:]+\
                  skq*(evk[1,3,:].conjugate()*evq[3,1,:]-evk[3,3,:].conjugate()*evq[1,1,:]))
    #Sqx
    so[s::4,2]=0.5*(evk[0,2,:].conjugate()*evq[2,0,:]+evk[2,2,:].conjugate()*evq[0,0,:]+\
               skq*(evk[1,2,:].conjugate()*evq[3,0,:]+evk[3,2,:].conjugate()*evq[1,0,:]))
    so[s+1::4,2]=0.5*(evk[0,3,:].conjugate()*evq[2,0,:]+evk[2,3,:].conjugate()*evq[0,0,:]+\
                 skq*(evk[1,3,:].conjugate()*evq[3,0,:]+evk[3,3,:].conjugate()*evq[1,0,:]))
    so[s+2::4,2]=0.5*(evk[0,2,:].conjugate()*evq[2,1,:]+evk[2,2,:].conjugate()*evq[0,1,:]+\
                 skq*(evk[1,2,:].conjugate()*evq[3,1,:]+evk[3,2,:].conjugate()*evq[1,1,:]))
    so[s+3::4,2]=0.5*(evk[0,3,:].conjugate()*evq[2,1,:]+evk[2,3,:].conjugate()*evq[0,1,:]+\
                 skq*(evk[1,3,:].conjugate()*evq[3,1,:]+evk[3,3,:].conjugate()*evq[1,1,:]))
    if (params['qx'],params['qy']) in [(0,0),(params['L']/2,params['L']/2)]:
        so[0,0]=0.5*np.sum(-evk[0,0,:].conjugate()*evq[0,0,:]+evk[2,0,:].conjugate()*evq[2,0,:]+\
                    skq*(-evk[1,0,:].conjugate()*evq[1,0,:]+evk[3,0,:].conjugate()*evq[3,0,:]))+\
                0.5*np.sum(-evk[0,1,:].conjugate()*evq[0,1,:]+evk[2,1,:].conjugate()*evq[2,1,:]+\
                    skq*(-evk[1,1,:].conjugate()*evq[1,1,:]+evk[3,1,:].conjugate()*evq[3,1,:]))
    return so
