import h5py
import matplotlib.pyplot as pl
from scipy.linalg import eigh
import scipy as sc
import stagflux as sf
import os
import re
import code

class InputFileError(Exception):
    def __init__(self,errstr):
        self.message=errstr

def CheckStat(filename,Nsamp=1):
    hfile=h5py.File(filename,'r')
    sw=0
    for r in hfile:
        for d in hfile[r]:
            sw+=1
    sdat=sw/Nsamp
    stats=sc.zeros(Nsamp)
    astat=0
    sample=0
    try:
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
    except Exception as err:
        print(err)
    return sdat

def GetEigSys(filename,Nsamp=1,channel=''):
    hfile=h5py.File(filename,'r')
    attr=GetAttr(filename)
    if channel=='':
        channel=attr['channel']
    sdat=CheckStat(filename,Nsamp)
    dat=hfile["/rank-1/data-0"]
    N=sc.shape(dat)[0]/2
    L=attr['L']
    q=[float(attr['qx'])/attr['L'], float(attr['qy'])/attr['L']]
    shift=[attr['phasex']/2.0,attr['phasey']/2.0]
    H=sc.zeros([Nsamp,N,N],complex)
    O=sc.zeros([Nsamp,N,N],complex)
    HO=sc.zeros([Nsamp,2*N,2*N])
    E=sc.zeros([Nsamp,N])
    V=sc.zeros([Nsamp,N,N],complex)
    ns=0
    sample=0
    for g in hfile:
        for d in hfile[g]:
            if sample<sc.shape(H)[0]:
                dat=hfile["/{0}/{1}".format(g,d)]
                HO[sample,:,:]+=dat
                H[sample,:,:]+=dat[0:N,0:2*N:2]+1j*dat[0:N,1:2*N:2]
                O[sample,:,:]+=dat[N:2*N,0:2*N:2]+1j*dat[N:2*N,1:2*N:2]
                ns+=1
                if ns==sdat:
                    H[sample,:,:]=0.5*(H[sample,:,:]+sc.conj(H[sample,:,:].T))/ns
                    O[sample,:,:]=0.5*(O[sample,:,:]+sc.conj(O[sample,:,:].T))/ns
                    HO[sample,:,:]/=ns
                    ns=0
                    sample+=1
    H,O=sf.fixfermisigns(attr['L'],attr['L'],shift,q,H,O,channel)
    print('{0} pair of (H,O) matrices loaded, now diagonalize'.format(sc.shape(H)[0]))
    for s in range(sc.shape(H)[0]):
        E[s,:],V[s,:,:]=eigh(sc.squeeze(H[s,:,:]),sc.squeeze(O[s,:,:]))
    print('diagonalization finished')
    return H,O,E,V

def GetSpinonOverlap(filename,Nsamp=1,channel='',O=None,r=None):
    attrs=GetAttr(filename)
    if channel=='':
        channel=attrs['channel']
    if channel=='long':
        raise ValueError('Longitudinal channel not yet implemented')
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    shift=[attrs['phasex']/2.0,attrs['phasey']/2.0]
    phi=attrs['phi']
    neel=attrs['neel']
    if O==None:
        H,O,E,V=GetEigSys(filename,Nsamp,channel)
    if r==None:
        X,Y=sc.meshgrid(range(-L/2,L/2),range(-L/2,L/2))
        r=sc.column_stack([X.flatten(),Y.flatten()])
    return sf.transspinonoverlap(O,L,L,q,shift,phi,neel,r)

def GetSqAmpl(filename,Nsamp=1,channel='',V=None,O=None,r=sc.zeros((1,2)),rp=sc.zeros((1,2))):
    """
    For the transverse channel:
    Calculates and returns Sq(sample,n,r)=<q,r|q,n><q,n|Sqp|GS>.
    Coordinates are in the order Sq(sample,n,r).
    For the longitudinal channel: input r has no effect.
    Calculates and return Sq(sample,n)=|<q,n|Sqz|GS>|^2.
    """
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
        return sf.sqwlongamp(V,O,L,L,q,shift,phi,neel)
    elif channel=='trans':
        return sf.sqwtransamp(V,O,L,L,q,shift,phi,neel,r,rp)
    pass

def GetSq(filename,Nsamp=1):
    attrs=GetAttr(filename)
    hfile=h5py.File(filename,'r')
    sdat=CheckStat(filename,Nsamp)
    filetype=''
    try:
        filetype=attrs['type']
    except KeyError:
        mo=re.match('.*/?[0-9]+-(StatSpinStruct)\.h5',filename)
        filetype=mo.groups()[0]
    if filetype!='StatSpinStruct':
        raise InputFileError('\"{0}\" is not a static structure factor file')
    N=pow(attrs['L'],2)
    Sq=sc.zeros((Nsamp,3,N))
    sample=0
    ns=0
    for r in hfile:
        for d in hfile[r]:
            Sq[sample,:,:]+=hfile[r][d][0:3,0::2]+1j*hfile[r][d][0:3,1::2]
            ns+=1
            if ns==sdat:
                ns=0
                sample+=1
    qx,qy=sc.meshgrid(range(attrs['L']),range(attrs['L']))
    return qx.flatten(),qy.flatten(),Sq

def PlotTransSpinon(filename,Nsamp=1,V=None,O=None,E=None,S=None,\
                    r=sc.zeros((1,2)),fig=None,width=0.1,\
                    shift=0,w=None):
    attrs=GetAttr(filename)
    channel='trans'
    try:
        channel=attrs['channel']
    except:
        pass
    if channel=='long':
        raise NotImplementedError('longitudinal channel not supported')
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    if E==None:
        H,O,E,V=GetEigSys(filename,Nsamp=Nsamp,channel=channel)
    if w==None:
        w=sc.arange(-6,6,0.01)
    if S==None:
        S=GetSqAmpl(filename,Nsamp=Nsamp,V=V,O=O,r=r)
    wqwr=sc.zeros((sc.shape(E)[0],sc.shape(w)[0]))
    for sample in range(Nsamp):
        En,Enp=sc.meshgrid(E[sample,:],E[sample,:])
        for ir in range(sc.shape(r)[0]):
            W,Wp=sc.meshgrid(sc.conj(S[sample,ir,:]),S[sample,ir,:])
            wqwr[sample,:]+=sf.gaussians(w,En.flatten()-Enp.flatten(),\
                                         W.flatten()*Wp.flatten(),\
                                         sc.ones(sc.shape(En.flatten()))*width)
    wqwr=wqwr/sc.shape(r)[0]+shift
    ax=None
    if fig is None:
        fig=pl.figure()
        ax=fig.gca()
    else:
        ax=fig.gca()
    for s in range(sc.shape(wqwr)[0]):
        ax.plot(w,wqwr[s,:])
    ax.set_xlim((sc.amin(w),sc.amax(w)))
    return fig

def GetAttr(filename):
    hfile=h5py.File(filename,'r')
    return hfile.attrs

def PlotSqw(filename,gsen,Nsamp=1,channel='',\
            fig=None,width=0.1,shift=0,\
            V=None,O=None,E=None,S=None,w=None):
    attrs=GetAttr(filename)
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    if E==None:
        H,O,E,V=GetEigSys(filename,Nsamp=Nsamp,channel=channel)
    if S==None:
        S=GetSqAmpl(filename,Nsamp=Nsamp,channel=channel,V=V,O=O)
    if w==None:
        w=sc.arange(-0.5,6,0.01)
    sqw=sc.zeros((sc.shape(E)[0],sc.shape(w)[0]),dtype=S.dtype)
    ax=None
    if len(sc.shape(S))==3:
        S=S[:,0,:]
    for s in range(sc.shape(sqw)[0]):
        sqw[s]=sf.gaussians(w,sc.squeeze(E[s,:])-gsen*L*L,sc.squeeze(S[s,:]),sc.ones(sc.shape(E)[1])*width)
    sqw+=shift
    if fig is None:
        fig=pl.figure()
        ax=fig.gca()
    else:
        ax=fig.gca()
    #ax.hold(True)
    for s in range(sc.shape(sqw)[0]):
        ax.plot(w,sqw[s,:])
        ax.plot(E[s,:]-gsen*L*L,sc.zeros(sc.shape(E)[1]),'o')
    if Nsamp!=1:
        ax.plot(w,sc.sum(sqw,0)/sc.shape(sqw)[0],'k--',linewidth=3)
    ax.set_xlim((sc.amin(w),sc.amax(w)))
    return fig

def ScanDir(folder='.',keys=[],return_dict=False):
    out={}
    for f in os.listdir(folder):
        if re.match(".*\.h5",f) is not None:
            out[f]=dict(GetAttr("{0}/{1}".format(folder,f)))
            s=f
            if len(keys):
                s="{0}: ".format(f)
                if keys=='*':
                    keys=out[f].keys()
                for k in keys:
                    try:
                        s="{0} {1} /".format(s,out[f][k])
                    except KeyError:
                        s="{0} None /".format(s)
            print(s)
    if return_dict:
        return out

if __name__=='__main__':
    PlotSqw('8-ProjHeis.h5',-0.664)
