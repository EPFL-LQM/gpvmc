import h5py
import matplotlib.pyplot as pl
from scipy.linalg import eigh
import scipy as sc
import stagflux as sf
import os
import re
import code

def GetEigSys(filename,Nsamp=1,channel=''):
    hfile=h5py.File(filename,'r')
    attr=GetAttr(filename)
    if channel=='':
        channel=attr['channel']
    sw=0
    for r in hfile:
        for d in hfile[r]:
            sw+=1
    sdat=sw/Nsamp
    stats=sc.zeros(Nsamp)
    astat=0
    sample=0
    for r in hfile:
        for d in hfile[r]:
            if sample<Nsamp:
                stats[sample]+=hfile[r][d].attrs['statistics']
                astat+=1
                if astat==sdat:
                    sample+=1
                    astat=0
    print("""{0} samples will be left out.
Average statistics of {1} (min: {2}, max: {3})"""\
            .format(sc.mod(sw,Nsamp),sc.mean(stats),\
            sc.amin(stats),sc.amax(stats)))
    dat=hfile["/rank-1/data-0"]
    N=sc.shape(dat)[0]/2
    L=attr['L']
    q=[float(attr['qx'])/attr['L'], float(attr['qy'])/attr['L']]
    shift=[attr['phasex']/2.0,attr['phasey']/2.0]
    H=sc.zeros([Nsamp,N,N],complex,'C')
    O=sc.zeros([Nsamp,N,N],complex,'C')
    E=sc.zeros([Nsamp,N])
    V=sc.zeros([Nsamp,N,N],complex)
    ns=0
    sample=0
    for g in hfile:
        for d in hfile[g]:
            if sample<sc.shape(H)[0]:
                dat=hfile["/{0}/{1}".format(g,d)]
                H[sample,:,:]+=dat[0:N,0:2*N:2]+1j*dat[0:N,1:2*N:2]
                O[sample,:,:]+=dat[N:2*N,0:2*N:2]+1j*dat[N:2*N,1:2*N:2]
                ns+=1
                if ns==sdat:
                    H[sample,:,:]=0.5*(H[sample,:,:]+sc.conj(H[sample,:,:].T))/ns
                    O[sample,:,:]=0.5*(O[sample,:,:]+sc.conj(O[sample,:,:].T))/ns
                    ns=0
                    sample+=1
    H,O=sf.fixfermisigns(attr['L'],attr['L'],shift,q,H,O,channel)
    for s in range(sc.shape(H)[0]):
        E[s,:],V[s,:,:]=eigh(sc.squeeze(H[s,:,:]),sc.squeeze(O[s,:,:]))
    return H,O,E,V

def GetSqAmpl(filename,Nsamp=1,channel='',V=None,O=None):
    attrs=GetAttr(filename)
    if channel=='':
        channel=attrs['channel']
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    shift=[attrs['phasex']/2.0,attrs['phasey']/2.0]
    phi=attrs['phi']
    neel=attrs['neel']
    if O==None or V== None:
        H,O,E,V=GetEigSys(filename,Nsamp,channel)
    if channel=='long':
        return E,sf.sqwlongamp(V,O,L,L,q,shift,phi,neel)
    elif channel=='trans':
        return E,sf.sqwtransamp(V,O,L,L,q,shift,phi,neel)
    pass

def GetTransSpinon(filename,time=sc.array([0]),Nsamp=1,channel='',V=None,O=None):
    attrs=GetAttr(filename)
    if channel=='':
        channel=attrs['channel']
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    shift=[attrs['phasex']/2.0,attrs['phasey']/2.0]
    phi=attrs['phi']
    neel=attrs['neel']
    if V==None or O==None:
        H,O,E,V=GetEigSys(filename,Nsamp,channel)
    X,Y=sc.meshgrid(sc.arange(-L/2,L/2+1,1),sc.arange(-L/2,L/2+1,1))
    r=sc.column_stack((X.flatten(),Y.flatten()))
    Wqrt=sf.transspinon(E,V,O,L,L,q,shift,phi,neel,r,time)
    return r, Wqrt

def GetAttr(filename):
    hfile=h5py.File(filename,'r')
    return hfile.attrs

def PlotSqw(filename,gsen,Nsamp=1,channel='',fig=None,width=0.1):
    attrs=GetAttr(filename)
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    E,S=GetSqAmpl(filename,Nsamp=Nsamp,channel=channel)
    w=sc.arange(-0.5,6,0.01)
    sqw=sc.zeros((sc.shape(E)[0],sc.shape(w)[0]))
    ax=None
    for s in range(sc.shape(sqw)[0]):
        sqw[s,:]=sf.gaussians(w,E[s,:]-gsen*L*L,S[s,:],sc.shape(E)[1]*[width])
    if fig is None:
        fig=pl.figure()
        ax=fig.gca()
    else:
        ax=fig.gca()
    print(ax)
    #ax.hold(True)
    for s in range(sc.shape(sqw)[0]):
        ax.plot(w,sqw[s,:])
        ax.plot(E[s,:]-gsen*L*L,sc.zeros(sc.shape(E)[1]),'o')
    return fig

def ScanDir(folder='.',keys=[]):
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
