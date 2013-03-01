import h5py
import matplotlib.pyplot as pl
import scipy as sc
import stagflux as sf
import os
import re

def GetEigSys(filename,channel=''):
    hfile=h5py.File(filename,'r')
    attr=GetAttr(filename)
    if channel=='':
        channel=attr['channel']
    dat=hfile["/rank-1/data-0"]
    N=sc.shape(dat)[0]/2
    L=attr['L']
    q=[float(attr['qx'])/attr['L'], float(attr['qy'])/attr['L']]
    shift=[attr['phasex']/2.0,attr['phasey']/2.0]
    H=sc.zeros([N,N],complex,'C')
    O=sc.zeros([N,N],complex,'C')
    ns=0
    for g in hfile:
        for d in hfile[g]:
            dat=hfile["/{0}/{1}".format(g,d)]
            H+=dat[0:N,0:2*N:2]+1j*dat[0:N,1:2*N:2]
            O+=dat[N:2*N,0:2*N:2]+1j*dat[N:2*N,1:2*N:2]
            ns+=1
    H=0.5/ns*(H+sc.conj(H.T))
    O=0.5/ns*(O+sc.conj(O.T))
    H,O=sf.fixfermisigns(attr['L'],attr['L'],shift,q,H,O,channel)
    return H,O

def GetSqAmpl(filename,channel=''):
    attrs=GetAttr(filename)
    if channel=='':
        channel=attrs['channel']
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    shift=[attrs['phasex']/2.0,attrs['phasey']/2.0]
    phi=attrs['phi']
    neel=attrs['neel']
    H,O=GetEigSys(filename,channel)
    if channel=='long':
        return sf.sqwlongamp(H,O,L,L,q,shift,phi,neel)
    elif channel=='trans':
        return sf.sqwtransamp(H,O,L,L,q,shift,phi,neel)
    pass

def GetTransSpinon(filename,time=sc.array([0]),channel=''):
    attrs=GetAttr(filename)
    if channel=='':
        channel=attrs['channel']
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    shift=[attrs['phasex']/2.0,attrs['phasey']/2.0]
    phi=attrs['phi']
    neel=attrs['neel']
    H,O=GetEigSys(filename,channel)
    X,Y=sc.meshgrid(sc.arange(-L/2,L/2+1,1),sc.arange(-L/2,L/2+1,1))
    r=sc.column_stack((X.flatten(),Y.flatten()))
    Wqrt=sf.transspinon(H,O,L,L,q,shift,phi,neel,r,time)
    nW=1.0/sc.sum(Wqrt,axis=0)
    return r, sc.einsum('ij,j->ij',Wqrt,nW)

def GetAttr(filename):
    hfile=h5py.File(filename,'r')
    return hfile.attrs

def PlotSqw(filename,gsen,channel='',fig=None):
    attrs=GetAttr(filename)
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    E,S=GetSqAmpl(filename,channel)
    w=sc.arange(-0.5,6,0.01)
    sqw=sf.gaussians(w,E-gsen*L*L,S,len(E)*[0.1])
    if fig is None:
        fig=pl.figure()
        ax=fig.gca()
    else:
        ax=fig.gca()
    ax.plot(w,sqw)
    ax.plot(E-gsen*L*L,sc.zeros(sc.shape(E)),'o')
    return fig

def ScanDir(folder,keys=[]):
    out={}
    for f in os.listdir(folder):
        if re.match(".*\.h5",f) is not None:
            print(f)
            out[f]=dict(GetAttr("{0}/{1}".format(folder,f)))
            if len(keys):
                s="{0}: ".format(f)
                for k in keys:
                    try:
                        s="{0} {1} /".format(s,out[f][k])
                    except KeyError:
                        s="{0} None /".format(s)
                print(s)
    return out

if __name__=='__main__':
    PlotSqw('8-ProjHeis.h5',-0.664)
