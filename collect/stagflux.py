import scipy as sc
import code

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

def phik(kx,ky,qx,qy,p):
    kqx=fold(kx-qx)
    kqy=fold(ky-qy)
    kqxm,kqym=mbzmod(kqx,kqy)
    #return -inmbz(kqx,kqy)*(  sc.conj(vk(kx,ky,1,-1,p))*vk(kqx,kqy,-1,1,p)\
    #                         +sc.conj(uk(kx,ky,1,-1,p))*uk(kqx,kqy,-1,1,p))\
    #       +(~inmbz(kqx,kqy))*(-sc.conj(vk(kx,ky,1,-1,p))*uk(kqx,kqy,-1,1,p)\
    #                           +sc.conj(uk(kx,ky,1,-1,p))*vk(kqx,kqy,-1,1,p))
    return  sc.conj(uk(kx,ky,1,1,p))*uk(kqx,kqy,-1,-1,p)+\
            sc.conj(vk(kx,ky,1,1,p))*vk(kqx,kqy,-1,-1,p)
    #return inmbz(kqx,kqy)*(sc.conj(vk(kx,ky,1,-1,p))*vk(kqxm,kqym,-1,1,p)\
    #                      +sc.conj(uk(kx,ky,1,-1,p))*uk(kqxm,kqym,-1,1,p))\
    #       +(~inmbz(kqx,kqy))*(sc.conj(vk(kx,ky,1,-1,p))*vk(kqxm,kqym,-1,1,p)\
    #                          -sc.conj(uk(kx,ky,1,-1,p))*uk(kqxm,kqym,-1,1,p))

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

def fermisigns(Lx,Ly,shift,q):
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

def sqwamp(H,O,Lx,Ly,q,shift,p):
    kx,ky=fermisea(Lx,Ly,shift)
    E,V=sc.linalg.eigh(H,O)
    sqw=len(E)*[0]
    for e in range(len(E)):
        pk=phik(kx,ky,q[0],q[1],[p['phi'],p['neel']])
        sqw[e]=abs(sc.dot(sc.conj(V[:,e]),sc.dot(O,pk)))**2
    return E,sqw

def gaussians(x,x0,A,sig):
    gg=sc.zeros(sc.shape(x))
    for p in range(len(x0)):
        gg+=A[p]*sc.sqrt(1/(2*sc.pi))/sig[p]*sc.exp(-0.5*(x-x0[p])**2/sig[p]**2)
    return gg

def Renorm(sqsq,O,Lx,Ly,q,shift,p):
    kx,ky=fermisea(Lx,Ly,shift)
    pp=phik(kx,ky,float(q[0])/Lx,float(q[1])/Ly,[p['phi'],p['neel']])
    b=sc.dot(sc.conj(pp),sc.dot(O,pp))
    r=sqsq/b
    if sc.isnan(r):
        r=1
    return r,b

if __name__=='__main__':
    pass
