import h5py
import matplotlib.pyplot as pl
import scipy as sc
import stagflux as sf

def GetEigSys(filename):
    hfile=h5py.File(filename,'r')
    N=int((hfile.attrs['L'])**2)
    H=sc.zeros([N/2,N/2],complex,'C')
    O=sc.zeros([N/2,N/2],complex,'C')
    ns=0
    for g in hfile:
        for d in hfile[g]:
            dat=hfile["/{0}/{1}".format(g,d)]
            H+=dat[0:N/2,0:N:2]+1j*dat[0:N/2,1:N:2]
            O+=dat[N/2:N,0:N:2]+1j*dat[N/2:N,1:N:2]
            ns+=1
    H=0.5/ns*(H+sc.conj(H.T))
    O=0.5/ns*(O+sc.conj(O.T))
    return H,O

def GetAttr(filename):
    hfile=h5py.File(filename,'r')
    return hfile.attrs

def PlotSqw(filename,gsen):
    attrs=GetAttr(filename)
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    shift=[attrs['phasex']/2.0,attrs['phasey']/2.0]
    H,O=GetEigSys(filename)
    [e,s]=sf.sqwamp(H,O,L,L,q,shift,attrs)
    w=sc.arange(-0.5,6,0.01)
    sqw=sf.gaussians(w,e-gsen*L*L,s,len(e)*[0.1])
    ax=pl.axes()
    ax.plot(w,sqw)
    ax.plot(e-gsen*L*L,sc.zeros(sc.shape(e)),'o')

if __name__=='__main__':
    PlotSqw('5-ProjHeis.h5',-0.664)
