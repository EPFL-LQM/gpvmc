import h5py
import matplotlib.pyplot as pl
import scipy as sc
import stagflux as sf

def GetEigSys(filename):
    hfile=h5py.File(filename,'r')
    attr=GetAttr(filename)
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
    H,O=sf.fixfermisigns(attr['L'],attr['L'],shift,q,H,O,attr['channel'])
    return H,O

def GetSqAmpl(filename):
    attrs=GetAttr(filename)
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    shift=[attrs['phasex']/2.0,attrs['phasey']/2.0]
    phi=attrs['phi']
    neel=attrs['neel']
    H,O=GetEigSys(filename)
    if attrs['channel']=='long':
        return sf.sqwlongamp(H,O,L,L,q,shift,phi,neel)
    elif attrs['channel']=='trans':
        return sf.sqwtransamp(H,O,L,L,q,shift,phi,neel)
    pass

def GetAttr(filename):
    hfile=h5py.File(filename,'r')
    return hfile.attrs

def PlotSqw(filename,gsen):
    attrs=GetAttr(filename)
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    E,S=GetSqAmpl(filename)
    w=sc.arange(-0.5,6,0.01)
    sqw=sf.gaussians(w,E-gsen*L*L,S,len(E)*[0.1])
    ax=pl.axes()
    ax.plot(w,sqw)
    ax.plot(E-gsen*L*L,sc.zeros(sc.shape(E)),'o')

if __name__=='__main__':
    PlotSqw('8-ProjHeis.h5',-0.664)
