import h5py
import warnings
import matplotlib.pyplot as pl
import matplotlib.mlab as ml
from scipy.linalg import eigh
import scipy as sc
from numpy.linalg import matrix_rank
import vmc_utils as vln
import stagflux as sf
import os
import re
import code

class InputFileError(Exception):
    def __init__(self,errstr):
        self.value=errstr
    def __str__(self):
        return str(self.value)

def GetStat(filename,Nsamp=1):
    hfile=[]
    if type(filename)==list:
        for i,f in enumerate(filename):
            hfile.append(h5py.File(f,'r'))
    elif type(filename)==str:
        hfile.append(h5py.File(filename,'r'))
        filename=[filename]
    stats=[]
    datapath=[]
    for ih,h in enumerate(hfile):
        for r in h:
            for d in h[r]:
                try:
                    stats.append(int(h[r][d].attrs['statistics'][0]))
                except KeyError as err:
                    stats.append(1)
                datapath.append((filename[ih],"/{0}/{1}".format(r,d)))
    bunches,args=vln.bunch(stats,Nsamp,indices=True)
    addstat=sc.array([sum(bunches[i]) for i in range(Nsamp)])
    print("Average statistics of {0} (min: {1}, max: {2})"\
            .format(sc.mean(addstat),sc.amin(addstat),sc.amax(addstat)))
    return datapath,args

def concatenate(infiles,outfile,Nsamp):
    if type(infiles)==str:
        infiles=[infiles]
    attr=GetAttr(infiles[0])
    dpath,args=GetStat(infiles,Nsamp)
    hout=h5py.File(outfile,'w')
    for k in attr:
        hout.attrs.create(k,attr[k])
    hfile=h5py.File(dpath[0][0],'r')
    dat=hfile['/rank-1/data-0']
    dtype=dat.dtype
    dshape=dat.shape
    for sample, bunch in enumerate(args):
        sdat=sc.zeros(dshape,dtype=dtype)
        st=0
        for d in bunch:
            hfile=h5py.File(dpath[d][0],'r')
            dat=hfile[dpath[d][1]]
            try:
                st+=dat.attrs['statistics']
            except:
                st+=1
            sdat+=dat
        sdat/=len(bunch)
        hdat=hout.create_dataset('/rank-1/data-{0}'.format(sample),dshape,dtype,sdat)
        hdat.attrs.create('statistics',st)
    hout.attrs.create('orig_files',infiles)
    hout.close()


def GetFermiSigns(filename,refstate=None,channel=None):
    attr=GetAttr(filename)
    filetype=''
    try:
        filetype=attr['type']
    except KeyError:
        mo=re.match('.*/?[0-9]+-(.*)\.h5',filename)
        filetype=mo.groups()[0]
    L=attr['L']
    if refstate==None:
        refstate=sc.zeros(2*L*L)
        refstate[0::2]=1
    if filetype=='WaveFunction':
        hfile=h5py.File(filename,'r')
        states=sc.column_stack([hfile['states_up'],hfile['states_do']])
        return sf.fermisigns(states,refstate)
    else:
        if channel==None:
            channel=attr['channel']
        L=attr['L']
        shift=[attr['phasex']/2.0,attr['phasey']/2.0]
        q=[float(attr['qx'])/L,float(attr['qy'])/L]
        if channel=='trans':
            return sf.transfermisigns(L,L,shift,q)
        elif channel=='long':
            return sf.longfermisigns(L,L,shift,q)
        else:
            raise KeyError('\"channel\" must be either \"trans\" or \"long\".')

def GetEigSys(filename,gsfile=None,Nsamp=1,channel=None,wavefile=None,q=None):
    if type(filename)==str:
        filename=[filename]
    hfile=h5py.File(filename[0],'r')
    attr=GetAttr(filename[0])
    if channel==None:
        channel=attr['channel']
    dpath,args=GetStat(filename,Nsamp)
    dat=hfile["/rank-1/data-0"]
    N=sc.shape(dat)[0]/2
    L=attr['L']
    shift=[attr['phasex']/2.0,attr['phasey']/2.0]
    H=sc.zeros([Nsamp,N,N],complex)
    O=sc.zeros([Nsamp,N,N],complex)
    E=sc.zeros([Nsamp,N])
    V=sc.zeros([Nsamp,N,N],complex)
    for sample,b in enumerate(args):
        for d in b:
            hfile=h5py.File(dpath[d][0],'r')
            dat=hfile[dpath[d][1]]
            H[sample,:,:]+=dat[0:N,0:2*N:2]+1j*dat[0:N,1:2*N:2]
            O[sample,:,:]+=dat[N:2*N,0:2*N:2]+1j*dat[N:2*N,1:2*N:2]
        H[sample,:,:]=0.5*(H[sample,:,:]+sc.conj(H[sample,:,:].T))/len(b)
        O[sample,:,:]=0.5*(O[sample,:,:]+sc.conj(O[sample,:,:].T))/len(b)
    if channel=='groundstate':
        return H
    fs=None
    refstate=sc.zeros(2*L*L)
    refstate[0::2]=1
    if wavefile==None:
        fs=GetFermiSigns(filename[0],refstate,channel=channel)
    else:
        fs=GetFermiSigns(wavefile,refstate,channel=channel)
    for s in range(sc.shape(H)[0]):
        H[s,:,:]=sc.dot(sc.diag(fs),sc.dot(H[s,:,:],sc.diag(fs)))
        O[s,:,:]=sc.dot(sc.diag(fs),sc.dot(O[s,:,:],sc.diag(fs)))
    ren=sc.ones(Nsamp)
    if gsfile!=None:
        ren=RenormalizeFactor(filename,gsfile,Nsamp=1,channel=channel,O=O,q=q)
    print('{0} pair of (H,O) matrices loaded, now diagonalize'.format(sc.shape(H)[0]))
    H=sc.einsum('ijk,i->ijk',H,ren)
    O=sc.einsum('ijk,i->ijk',O,ren)
    for s in range(sc.shape(H)[0]):
        E[s,:],V[s,:,:]=vln.geneigh(sc.squeeze(H[s,:,:]),sc.squeeze(O[s,:,:]))
    print('diagonalization finished')
    return H,O,E,V

def RenormalizeFactor(excfile,gsfile,channel=None,Nsamp=1,O=None,q=None):
    if type(excfile)==str:
        excfile=[excfile]
    if type(gsfile)==str:
        gsfile=[gsfile]
    exat=GetAttr(excfile[0])
    gsat=GetAttr(gsfile[0])
    L=exat['L']
    if q==None:
        q=sc.array([exat['qx'],exat['qy']])
    shift=sc.array([exat['phasex']/2.0,exat['phasey']/2.0])
    phi=exat['phi']
    neel=exat['neel']
    qx,qy,Sq=GetSq(gsfile)
    kx,ky=sf.fermisea(L,L,shift)
    qidx=ml.find((qx==q[0])*(qy==q[1]))
    if O==None:
        _,O,_,_=GetEigSys(excfile,Nsamp)
    pk=None
    sqq=None
    if channel==None:
        channel=exat['channel']
    if channel=='trans':
        pk=sc.squeeze(sf.phiktrans(kx,ky,q[0]/L,q[1]/L,[phi,neel]))
        sqq=sc.real(0.5*(Sq[0,1,qidx]+Sq[0,2,qidx]))
    elif channel=='long':
        pkup=sc.squeeze(sf.phiklong(kx,ky,q[0]/L,q[1]/L,1,[phi,neel]))
        pkdo=sc.squeeze(sf.phiklong(kx,ky,q[0]/L,q[1]/L,-1,[phi,neel]))
        if (q[0]/L==0.5 and q[1]/L==0.5) or (q[0]/L==0 and q[1]/L==0):
            pk=sc.zeros(2*sc.shape(pkup)[0]+1,complex)
        else:
            pk=sc.zeros(2*sc.shape(pkup)[0],complex)
        pk[0:2*sc.shape(pkup)[0]:2]=pkup
        pk[1:2*sc.shape(pkdo)[0]:2]=pkdo
        if (qx[0]/L==0.5 and q[1]/L==0.5) or (q[0]/L==0 and q[1]/L==0):
            if neel==0:
                pk[-1]=0
            else:
                pk[-1]=sum(neel/sf.omega(kx,ky,[phi,neel]))
        sqq=Sq[0,0,qidx]
    else:
        raise(InputFileError('In file \''+excfile+'\', channel=\''+str(channel)+'\'. Should be \'trans\' or \'long\''))
    sqe=sc.einsum('i,jik,k->j',sc.conj(pk),O,pk)
    out=sc.zeros(Nsamp)
    for n in range(Nsamp):
        if abs(sqq)<1e-6 or abs(sqe[n])<1e-6:
            warnings.warn('Probably ill-defined renormalization, returns 1 for sample {0} out of {1}'.format(n,Nsamp),UserWarning)
            out[n]=1
        else:
            out[n]=sc.real(sqq/sqe[n])
    return out

def GetSpinonOverlap(filename,Nsamp=1,channel=None,O=None,V=None,r=None):
    if type(filename)==str:
        filename=[filename]
    attrs=GetAttr(filename[0])
    if channel==None:
        channel=attrs['channel']
    if channel=='long':
        raise ValueError('Longitudinal channel not yet implemented')
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    shift=[attrs['phasex']/2.0,attrs['phasey']/2.0]
    phi=attrs['phi']
    neel=attrs['neel']
    if O==None:
        H,O,E,V=GetEigSys(filename,Nsamp=Nsamp,channel=channel)
    if r==None:
        X,Y=sc.meshgrid(range(-L/2,L/2),range(-L/2,L/2))
        r=sc.column_stack([X.flatten(),Y.flatten()])
    return sf.transspinonoverlap(O,V,L,L,q,shift,phi,neel,r)

def GetSqAmpl(filename,Nsamp=1,channel=None,V=None,O=None,r=sc.zeros((1,2)),rp=sc.zeros((1,2)),q=None):
    """
    For the transverse channel:
    Calculates and returns Sq(sample,n,r)=<q,r|q,n><q,n|Sqp|GS>.
    Coordinates are in the order Sq(sample,n,r).
    For the longitudinal channel: input r has no effect.
    Calculates and return Sq(sample,n)=|<q,n|Sqz|GS>|^2.
    """
    if type(filename)==str:
        filename=[filename]
    attrs=GetAttr(filename[0])
    if channel==None:
        channel=attrs['channel']
    L=attrs['L']
    if q==None:
        q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    else:
        q=[float(q[0])/L,float(q[1])/L]
    shift=[attrs['phasex']/2.0,attrs['phasey']/2.0]
    phi=attrs['phi']
    neel=attrs['neel']
    if O==None or V== None:
        H,O,E,V=GetEigSys(filename,Nsamp=Nsamp,channel=channel,q=q)
    if channel=='long':
        return sf.sqwlongamp(V,O,L,L,q,shift,phi,neel)
    elif channel=='trans':
        return sf.sqwtransamp(V,O,L,L,q,shift,phi,neel,r,rp)
    pass

def GetSq(filename,Nsamp=1):
    if type(filename)==str:
        filename=[filename]
    attrs=GetAttr(filename[0])
    #hfile=h5py.File(filename,'r')
    dpath,args=GetStat(filename,Nsamp)
    filetype=''
    try:
        filetype=attrs['type']
    except KeyError:
        mo=re.match('.*/?[0-9]+-(StatSpinStruct)\.h5',filename)
        filetype=mo.groups()[0]
    if filetype!='StatSpinStruct':
        raise InputFileError('\"{0}\" is not a static structure factor file')
    N=pow(attrs['L'],2)
    Sq=sc.zeros((Nsamp,3,N),complex)
    for sample,b in enumerate(args):
        for d in b:
            hfile=h5py.File(dpath[d][0],'r')
            Sq[sample,:,:]+=hfile[dpath[d][1]][0:3,0::2]+1j*hfile[dpath[d][1]][0:3,1::2]
        Sq[sample,:,:]/=len(b)
    qx,qy=sc.meshgrid(range(attrs['L']),range(attrs['L']))
    return qx.flatten(),qy.flatten(),Sq

def GetAttr(filename):
    hfile=h5py.File(filename,'r')
    return hfile.attrs

def PlotSqw(filename,gsen,Nsamp=1,channel=None,\
            fig=None,width=0.1,shift=0,\
            V=None,O=None,E=None,S=None,w=None,
            gsspinfile=None, wavefile=None):
    if type(filename)==str:
        filename=[filename]
    attrs=GetAttr(filename[0])
    L=attrs['L']
    q=[float(attrs['qx']/L),float(attrs['qy'])/L]
    if E==None:
        H,O,E,V=GetEigSys(filename,Nsamp=Nsamp,channel=channel,gsfile=gsspinfile,wavefile=wavefile)
    if S==None:
        S=GetSqAmpl(filename,Nsamp=Nsamp,channel=channel,V=V,O=O)
    if w==None:
        w=sc.arange(-0.5,6,0.01)
    sqw=sc.zeros((sc.shape(E)[0],sc.shape(w)[0]),dtype=S.dtype)
    ax=None
    S=S[:,0,0,:]
    for s in range(sc.shape(sqw)[0]):
        idx=~sc.isnan(E[s,:])
        sqw[s]=sf.gaussians(w,sc.squeeze(E[s,idx])-gsen*L*L,sc.squeeze(S[s,idx]),sc.ones(sc.shape(E[s,idx]))*width)
    sqw+=shift
    if fig is None:
        fig=pl.figure()
        ax=fig.gca()
    else:
        ax=fig.gca()
    #ax.hold(True)
    for s in range(sc.shape(sqw)[0]):
        ax.plot(w,sc.real(sqw[s,:]))
        ax.plot(E[s,:]-gsen*L*L,sc.squeeze(S[s,:])*sc.sqrt(1/2.0/sc.pi)/width,'o')
    if Nsamp!=1:
        ax.plot(w,sc.sum(sqw,0)/sc.shape(sqw)[0],'k--',linewidth=3)
    ax.set_xlim((sc.amin(w),sc.amax(w)))
    return fig

def ScanDir(folder='.',keys=[],pattern=r".*\.h5",return_dict=False):
    out={}
    for f in os.listdir(folder):
        if re.match(pattern,f) is not None:
            try:
                out[folder+'/'+f]=dict(GetAttr("{0}/{1}".format(folder,f)))
                s=f
                if len(keys):
                    s="{0}: ".format(f)
                    if keys=='*':
                        keys=out[folder+'/'+f].keys()
                    for k in keys:
                        try:
                            s="{0} {1}:{2} /".format(s,k,out[folder+'/'+f][k])
                        except KeyError:
                            s="{0} None /".format(s)
                print(s)
            except IOError:
                print('Could not open \"'+f+'\".')
    if return_dict:
        return out

if __name__=='__main__':
    PlotSqw('8-ProjHeis.h5',-0.664)
