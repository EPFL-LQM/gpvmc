import numpy as np

def mbz(params):
    Lx=int(params['L'])
    Ly=int(params['L'])
    shift=[params['phase_shift_x']/2.0,params['phase_shift_y']/2.0]
    fskx=np.zeros(Lx*Ly/2)
    fsky=np.zeros(Lx*Ly/2)
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
    return k-((k<0) + (k>=1))*np.floor(k)

def mbzmod(kx,ky):
    mkx=fold(kx)
    mky=fold(ky)
    mmkx=fold(mkx+0.5*(1-inmbz(mkx,mky)))
    mmky=fold(mky+0.5*(1-inmbz(mkx,mky)))
    return mmkx,mmky
