#!/bin/env python
import h5py
import scipy as sc
from scipy.linalg import eigh
import stagflux as sf

def check_sampling(hfile,Nsamp):
    samptot=0
    for r in hfile:
        samptot+=len(hfile[r])
    print('Will drop {0} out of {1} samples'.format(sc.mod(samptot,Nsamp),samptot))
    return samptot-sc.mod(samptot,Nsamp)

def GetEigSys(rootdir,L,q,calc_id,Nsamp):
    infile="{0}/L{1}x{1}/q_{2}_{3}/{4}-SpinExcitonEnergy.h5".format(\
            rootdir,L,q[0],q[1],calc_id)
    try:
        hfile=h5py.File(infile,'r')
    except IOError as err:
        print("could not open file \"{0}\".".format(infile))
        raise err
    samptot=check_sampling(hfile,Nsamp)
    N=int((hfile.attrs['L'])**2)
    H=sc.zeros([Nsamp,N/2,N/2],complex,'C')
    O=sc.zeros([Nsamp,N/2,N/2],complex,'C')
    ns=0
    i=0
    for g in hfile:
        for d in hfile[g]:
            if i<Nsamp:
                dat=hfile["/{0}/{1}".format(g,d)]
                H[i,:,:]+=dat[0:N/2,0:N:2]+dat[0:N/2,1:N:2]*1j
                O[i,:,:]+=dat[N/2:N,0:N:2]+dat[N/2:N,1:N:2]*1j
                ns+=1
                if ns==samptot/Nsamp:
                    H[i,:,:]/=ns
                    O[i,:,:]/=ns
                    i+=1
                    ns=0
    for i in range(Nsamp):
        H[i,:,:]=0.5*(H[i,:,:]+sc.conj(H[i,:,:].T))
        O[i,:,:]=0.5*(O[i,:,:]+sc.conj(O[i,:,:].T))
    p=GetParams(rootdir,L,q,calc_id)
    shift=[p['phasex']/2.0,p['phasey']/2.0]
    fs=sf.fermisigns(L,L,shift,q)
    H=sc.dot(sc.diag(fs),sc.dot(H,sc.diag(fs)))
    O=sc.dot(sc.diag(fs),sc.dot(O,sc.diag(fs)))
    return H, O

def GetParams(rootdir,L,q,calc_id):
    infile="{0}/L{1}x{1}/q_{2}_{3}/{4}-SpinExcitonEnergy.h5".format(\
            rootdir,L,q[0],q[1],calc_id)
    hfile=h5py.File(infile,'r')
    odic=dict()
    for a in hfile.attrs:
        odic[a]=hfile.attrs[a]
    return odic

def GetGroundState(rootdir,L,calc_id,Nsamp):
    infile="{0}/L{1}x{1}/GS/{2}-HeisEn.h5".format(\
            rootdir,L,calc_id)
    hfile=h5py.File(infile,'r')
    samptot=check_sampling(hfile,Nsamp)
    GSen=sc.zeros([Nsamp],complex,'C')
    ns=0
    i=0
    for g in hfile:
        for d in hfile[g]:
            dat=hfile["{0}/{1}".format(g,d)]
            GSen[i]+=dat[0,0]+dat[0,1]*1j
            ns+=1
            if ns==samptot/Nsamp:
                GSen[i]/=ns
                i+=1
                ns=0
    return GSen

def GetSingleMode(rootdir,L,calc_id,Nsamp):
    infile="{0}/L{1}x{1}/GS/{2}-SpinDensity.h5".format(\
            rootdir,L,calc_id)
    hfile=h5py.File(infile,'r')
    samptot=check_sampling(hfile,Nsamp)
    Nqs=sc.shape(hfile['rank-1/data-0'])[1]/2
    Sq=sc.zeros([Nsamp,Nqs],complex,'C')
    Hq=sc.zeros([Nsamp,Nqs],complex,'C')
    ns=0
    i=0
    for g in hfile:
        for d in hfile[g]:
            dat=hfile["{0}/{1}".format(g,d)]
            Hq[i,:]+=dat[0,0::2]+dat[0,1::2]*1j
            Sq[i,:]+=dat[1,0::2]+dat[1,1::2]*1j
            ns+=1
            if ns==samptot/Nsamp:
                Sq[i]/=ns
                Hq[i]/=ns
                i+=1
                ns=0
    return Hq,Sq

if __name__=='__main__':
    [H,O]=GetEigSys(rootdir='/home/bastien/data/VarMC/local/testing',\
            L=8,q=[2,2],calc_id=40,Nsamp=1)
