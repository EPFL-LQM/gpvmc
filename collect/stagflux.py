import scipy as sc
from scipy.linalg import eigh
import code
import warnings

def delta(kx,ky,param):
    return 0.5*(sc.cos(kx*2*sc.pi)*sc.exp(1j*param[0]*sc.pi)\
            +sc.cos(ky*2*sc.pi)*sc.exp(-1j*param[0]*sc.pi))

def omega(kx,ky,param):
    return sc.sqrt(param[1]**2+abs(delta(kx,ky,param))**2)

def uk(kx,ky,spin,band,param):
    if band==-1:
        return sc.sqrt(0.5*(1+spin*param[1]/omega(kx,ky,param)))
    else:
        return  -sc.conj(delta(kx,ky,param))/abs(delta(kx,ky,param))\
                *sc.sqrt(0.5*(1-spin*param[1]/omega(kx,ky,param)))

def vk(kx,ky,spin,band,param):
    if band==-1:
        return  delta(kx,ky,param)/abs(delta(kx,ky,param))*\
                sc.sqrt(0.5*(1-spin*param[1]/omega(kx,ky,param)));
    else:
        return sc.sqrt(0.5*(1+spin*param[1]/omega(kx,ky,param)));

def phiktrans(kx,ky,qx,qy,p,r=sc.zeros((1,2))):
    """
    Returns phi[k,r] such that |q,r>=sum_k phi[k,r]|q,k>
    """
    kqx=kx-qx
    kqy=ky-qy
    pk=sc.zeros((sc.shape(kx)[0],sc.shape(r)[0]),complex)
    pke=sc.conj(uk(kx,ky,1,1,p))*uk(kqx,kqy,-1,-1,p)+\
        sc.conj(vk(kx,ky,1,1,p))*vk(kqx,kqy,-1,-1,p)
    pko=sc.conj(vk(kx,ky,1,1,p))*uk(kqx,kqy,-1,-1,p)+\
        sc.conj(uk(kx,ky,1,1,p))*vk(kqx,kqy,-1,-1,p)
    even=1-sc.mod(r[:,0]+r[:,1],2)
    odd=sc.mod(r[:,0]+r[:,1],2)
    ph=sc.exp(-2j*sc.pi*(sc.einsum('i,j->ij',kx,r[:,0])+sc.einsum('i,j->ij',ky,r[:,1])))
    pk=sc.einsum('ij,j,i->ij',ph,even,pke)+sc.einsum('ij,j,i->ij',ph,odd,pko)
    return pk

def phiklong(kx,ky,qx,qy,spin,p):
    kqx=kx-qx
    kqy=ky-qy
    pk=0.5*spin*(sc.conj(uk(kx,ky,spin,1,p))*uk(kqx,kqy,spin,-1,p)+\
                 sc.conj(vk(kx,ky,spin,1,p))*vk(kqx,kqy,spin,-1,p))
    return pk

def fermisea(Lx,Ly,shift):
    fskx=sc.zeros(Lx*Ly/2)
    fsky=sc.zeros(Lx*Ly/2)
    i=0
    for kx in range(Lx):
        for ky in range(Ly):
            if inmbz(float(kx+shift[0])/Lx,float(ky+shift[1])/Ly):
                fskx[i]=float(kx+shift[0])/Lx
                fsky[i]=float(ky+shift[1])/Ly
                i+=1
    return fskx,fsky

def inmbz(kx,ky):
    fkx=fold(kx)
    fky=fold(ky)
    tol=1e-6
    return (fky<=-0.5+fkx+tol) + (fky< 0.5-fkx-tol)\
            + (fky>=1.5-fkx-tol) + (fky>0.5+fkx+tol)

def fold(k):
    return k-((k<0) + (k>=1))*sc.floor(k)

def mbzmod(kx,ky):
    mkx=fold(kx)
    mky=fold(ky)
    mmkx=fold(mkx+0.5*(1-inmbz(mkx,mky)))
    mmky=fold(mky+0.5*(1-inmbz(mkx,mky)))
    return mmkx,mmky

def shiftmbz(kx,ky):
    kx,ky=mbzmod(kx,ky)
    mkx=sc.array([k-(k>0.5) for k in kx])
    mky=sc.array([k-(k>0.5) for k in ky])
    return mkx,mky

def transfermisigns(Lx,Ly,shift,q):
    kx,ky=fermisea(Lx,Ly,shift)
    kqx=sc.mod(kx*Lx-shift[0]-q[0]*Lx,Lx)/Lx+shift[0]/Lx
    kqy=sc.mod(ky*Ly-shift[1]-q[1]*Ly,Ly)/Ly+shift[1]/Ly
    kqx,kqy=mbzmod(kqx,kqy)
    fsign=sc.zeros(len(kx))
    gsup=sc.zeros(2*len(kx))
    gsup[0::2]=1
    gsdo=sc.zeros(2*len(kx))
    gsdo[0::2]=1
    for k in range(len(kx)):
        ma=abs(kqx[k]-kx)+abs(kqy[k]-ky)
        ma=(ma==sc.amin(ma))
        idx=sc.arange(0,len(ma),1)[ma]
        fsign[k]=(-1)**sum(gsdo[0:2*idx])*(-1)**sum(gsup[0:(2*k+1)])
    return fsign

def longfermisigns(Lx,Ly,shift,q):
    kx,ky=fermisea(Lx,Ly,shift)
    kqx=sc.mod(kx*Lx-shift[0]-q[0]*Lx,Lx)/Lx+shift[0]/Lx
    kqy=sc.mod(ky*Ly-shift[1]-q[1]*Ly,Ly)/Ly+shift[1]/Ly
    kqx,kqy=mbzmod(kqx,kqy)
    fsign=sc.zeros(2*len(kx))
    if (abs(q[0])+abs(q[1]))<1e-6 or (abs(q[0]-0.5)+abs(q[1]-0.5))<1e-6:
        fsign=sc.zeros(2*len(kx)+1)
        fsign[-1]=1
    for k in range(len(kx)):
        gsup=sc.zeros(2*len(kx))
        gsup[0::2]=1
        ma=abs(kqx[k]-kx)+abs(kqy[k]-ky)
        ma=(ma==sc.amin(ma))
        idx=sc.arange(0,len(ma),1)[ma]
        fsign[2*k]=(-1)**sum(gsup[0:2*idx])
        gsup[2*idx]=0
        fsign[2*k]*=(-1)**sum(gsup[0:(2*k+1)])
        fsign[2*k+1]=fsign[2*k]
    return fsign

def fermisigns(states,refstate):
    fs=sc.ones(sc.shape(states)[0])
    for s in range(sc.shape(states)[0]):
        hole=None
        part=None
        ref=sc.array(refstate)
        for i in range(sc.shape(states)[1]):
            if ref[i]==0 and states[s,i]!=0:
                part=i
            if ref[i]!=0 and states[s,i]==0:
                hole=i
            if part!=None and hole!=None:
                fs[s]*=(-1)**(sum(ref[min(hole,part)+1:max(hole,part)]))
                ref[hole]=0
                ref[part]=1
                hole=None
                part=None
    return fs

def fixfermisigns(Lx,Ly,shift,q,H,O,ori):
    fs=[]
    if ori=='trans':
        fs=transfermisigns(Lx,Ly,shift,q)
    elif ori=='long':
        fs=longfermisigns(Lx,Ly,shift,q)
    H=sc.einsum('i,jik,k->jik',fs,H,fs)
    O=sc.einsum('i,jik,k->jik',fs,O,fs)
    return H,O

def sqwtransamp(V,O,Lx,Ly,q,shift,phi,neel,r=sc.zeros((1,2)),rp=sc.zeros((1,2))):
    """
    Returns Sq[sample,r,n]=<q,r|q,n><q,n|q,0>
    """
    sqn=sc.zeros(sc.shape(V)[0:2],complex)
    kx,ky=fermisea(Lx,Ly,shift)
    pkrp=sc.zeros((sc.shape(V)[1],sc.shape(rp)[0]))
    pkr=sc.zeros((sc.shape(V)[1],sc.shape(rp)[0]))
    pkrp[0:len(kx),:]=phiktrans(kx,ky,q[0],q[1],[phi,neel],rp)
    pkr[0:len(ky),:]=phiktrans(kx,ky,q[0],q[1],[phi,neel],r)
    OV=sc.einsum('ijk,ikl->ijl',O,V)
    rhs=sc.einsum('ijk,jl->ikl',sc.conj(OV),pkrp)
    lhs=sc.einsum('ij,kil->kjl',sc.conj(pkr),OV)
    sqn=sc.einsum('ijk,ikj->ijk',lhs,rhs)
    return sqn

def sqwlongamp(V,O,Lx,Ly,q,shift,phi,neel):
    sqn=sc.zeros(sc.shape(V)[0:2],complex)
    kx,ky=fermisea(Lx,Ly,shift)
    pkup=phiklong(kx,ky,q[0],q[1],1,[phi,neel])
    pkdo=phiklong(kx,ky,q[0],q[1],-1,[phi,neel])
    pk=sc.zeros(sc.shape(V)[1],complex)
    pk[0:2*len(pkup):2]=pkup
    pk[1:2*len(pkup):2]=pkdo
    if (abs(q[0])+abs(q[1]))<1e-6 or\
       (abs(q[0]-0.5)+abs(q[1]-0.5))<1e-6:
        pk=sc.append(pk,0)
        if neel!=0:
            pk[-1]=sc.sum(neel/omega(kx,ky,[phi,neel]))
    sqn=abs(sc.einsum('ijk,ijl,l->ik',sc.conj(V),O,pk))**2
    return sqn

def transspinonoverlap(O,V,Lx,Ly,q,shift,phi,neel,r):
    kx,ky=fermisea(Lx,Ly,shift)
    pkr=phiktrans(kx,ky,q[0],q[1],[phi,neel],r)
    ork=sc.einsum('ij,kil->kjl',sc.conj(pkr),O)
    return sc.einsum('kil,klm->kim',ork,V)

def gaussians(x,x0,A,sig):
    #if sc.amax(abs(sc.imag(A)))/sc.amax(abs(sc.real(A)))>0.01:
    #    warnings.warn(\
    #'Gaussian amplitude has a sizable imaginary part\(max(|Im|)/max(|Re|)={0}, mean(abs(A))={1}).'\
    #        .format(sc.amax(abs(sc.imag(A)))/sc.amax(abs(sc.real(A))), sc.mean(abs(A))))
    amp=A*sc.sqrt(1/2.0/sc.pi)/sig
    [X,X0]=sc.meshgrid(x,x0)
    gg=sc.einsum('i,ij',amp,sc.exp(-0.5*(X-X0)**2/sc.tile(sig**2,(sc.shape(x)[0],1)).T))
    return gg

def Renorm(sqsq,O,Lx,Ly,q,shift,p):
    kx,ky=fermisea(Lx,Ly,shift)
    pp=phiktrans(kx,ky,float(q[0])/Lx,float(q[1])/Ly,[p['phi'],p['neel']])
    b=sc.dot(sc.conj(pp),sc.dot(O,pp))
    r=sqsq/b
    if sc.isnan(r):
        r=1
    return r,b

if __name__=='__main__':
    pass
