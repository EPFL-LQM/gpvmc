#include "SFpNxpHzWaveFunction.h"
#include <cstdlib>
#include <stdexcept>
#include <iostream>

using namespace std;

SFpNxpHzWaveFunction::SFpNxpHzWaveFunction(size_t Lx, size_t Ly,
                                           size_t Nby,
                                           double phi, double neel, double hz,
                                           const vector<double>& bc_phase)
    :m_Lx(Lx), m_Ly(Ly),
     m_qn2fock(new size_t[Lx*Ly*4]),
     m_fock2qn(new size_t[Lx*Ly*2*3]),
     m_phi(phi), m_nx(neel), m_hz(hz),
     m_bc_phase(bc_phase)
{
    vector<size_t> Nbyv(1,Nby);
    vector<size_t> Nfs(1,2*Lx*Ly);
    build_base(Nbyv,Nfs);
    size_t idx=0;
    for(size_t i=0;i<2*Lx*Ly;++i) m_qn2fock[i]=Lx*Ly;
    for(size_t kx=0;kx<Lx;++kx){
        for(size_t ky=0;ky<Ly;++ky){
            if(inmbz(kx,ky)){
                for(size_t b=0;b<4;++b){
                    m_fock2qn[idx*3]=kx;
                    m_fock2qn[idx*3+1]=ky;
                    m_fock2qn[idx*3+2]=b;
                    m_qn2fock[(kx*Ly+ky)*4+b]=idx;
                    idx+=1;
                }
            }
        }
    }
    init_matrices();
}

SFpNxpHzWaveFunction::~SFpNxpHzWaveFunction()
{
    delete [] m_fock2qn;
    delete [] m_qn2fock;
}

std::complex<double> SFpNxpHzWaveFunction::matrix_element(size_t fk,size_t fr, size_t fl)
{
    std::complex<double> I(0.,1.);
    double k[2], r[2];
    k[0]=(m_fock2qn[3*fk]+m_bc_phase[0]/(2.0))*2*M_PI/double(m_Lx);
    k[1]=(m_fock2qn[3*fk+1]+m_bc_phase[1]/(2.0))*2*M_PI/double(m_Ly);
    size_t band=m_fock2qn[3*fk+2];
    int ifr=fr;
    size_t up=ifr%2;
    r[1]=((ifr-up)/2)%m_Ly;
    r[0]=((ifr-up)/2-r[1])/m_Ly;
    up=size_t(1-int(up));
    std::complex<double> out;
    double phir=1-int(r[0]+r[1])%2;
    if(up){
        out=exp(I*(k[0]*r[0]+k[1]*r[1]))*
            (phir*(Uk(k,band)+Xk(k,band))
             +(1-phir)*(Vk(k,band)+Yk(k,band)));
    } else {
        out=exp(I*(k[0]*r[0]+k[1]*r[1]))*
            (phir*(-Uk(k,band)+Xk(k,band))
             +(1-phir)*(-Vk(k,band)+Yk(k,band)));
    }
#ifdef CRAY
    if(isnan(real(out)) || isnan(imag(out)))
#else
    if(std::isnan(real(out)) || std::isnan(imag(out)))
#endif
    {
#ifdef EXCEPT
        throw std::logic_error("SFpNxpHzWaveFunction::matrix_element::NaN encoutered");
#else
        std::cerr<<"SFpNxpHzWaveFunction::matrix_element::NaN encoutered"<<std::endl;
        abort();
#endif
    }
    return out;
}

bool SFpNxpHzWaveFunction::inmbz(size_t kx, size_t ky) const
{
    double kkx=(double(kx)+m_bc_phase[0]/2.0)/m_Lx;
    double kky=(double(ky)+m_bc_phase[1]/2.0)/m_Ly;
    double eps=1e-9;
    bool out= (kkx>=0 && kkx<1 &&
            kky>=0 && kky<1) &&
           (kkx+kky-0.5<0 ||
            kkx+kky-1.5>=-eps ||
            kky-kkx-0.5>eps ||
            kky-kkx+0.5<=eps);
    return out;
}

void SFpNxpHzWaveFunction::mbzmod(size_t* k) const
{
    if(inmbz(k[0],k[1])) return;
    else if(k[0]<m_Lx/2 && k[1]<m_Ly/2){
        k[0]+=m_Lx/2;
        k[1]+=m_Ly/2;
    } else if(k[0]<m_Lx/2 && k[1]>=m_Ly/2){
        k[0]+=m_Lx/2;
        k[1]-=m_Ly/2;
    } else if(k[0]>=m_Lx/2 && k[1]<m_Ly/2){
        k[0]-=m_Lx/2;
        k[1]+=m_Ly/2;
    } else {
        k[0]-=m_Lx/2;
        k[1]-=m_Ly/2;
    }
}

int SFpNxpHzWaveFunction::sign(double x) const
{
    if(x>0)
        return 1;
    else if(x<0)
        return -1;
    else
        return 0;
}

std::complex<double> SFpNxpHzWaveFunction::delta(double* k) const
{
    return 0.5*(cos(k[0])*std::polar(1.0,m_phi)+cos(k[1])*std::polar(1.0,-m_phi));
}

double SFpNxpHzWaveFunction::omega(double* k, size_t band) const
{
    if(band==0){
        return -sqrt(m_nx*m_nx + pow(m_hz-abs(delta(k)),2));
    } else if(band==1){
        return -sqrt(m_nx*m_nx + pow(m_hz+abs(delta(k)),2));
    } else if(band==2){
        return sqrt(m_nx*m_nx + pow(m_hz-abs(delta(k)),2));
    } else {//band==3
        return sqrt(m_nx*m_nx + pow(m_hz+abs(delta(k)),2));
    }
}

std::complex<double> SFpNxpHzWaveFunction::Uk(double* k, size_t band) const
{
    double om=omega(k,band);
    return 0.5*sqrt(1.0-m_nx/om);
}

std::complex<double> SFpNxpHzWaveFunction::Vk(double* k, size_t band) const
{
    double om=omega(k,band);
    std::complex<double> d=delta(k);
    std::complex<double> pd=d/abs(d);
    if(band==0)
        return 0.5*sign(abs(d)-m_hz)*pd*sqrt(1+m_nx/om);
    if(band==1)
        return 0.5*sign(abs(d)+m_hz)*pd*sqrt(1+m_nx/om);
    if(band==2)
        return -0.5*sign(abs(d)-m_hz)*pd*sqrt(1+m_nx/om);
    else //band==3
        return -0.5*sign(abs(d)+m_hz)*pd*sqrt(1+m_nx/om);
}

std::complex<double> SFpNxpHzWaveFunction::Xk(double* k, size_t band) const
{
    std::complex<double> d=delta(k);
    double om=omega(k,band);
    if(band==0)
        return 0.5*sign(abs(d)-m_hz)*sqrt(1+m_nx/om);
    if(band==1)
        return -0.5*sign(abs(d)+m_hz)*sqrt(1+m_nx/om);
    if(band==2)
        return -0.5*sign(abs(d)-m_hz)*sqrt(1+m_nx/om);
    else //band==3
        return 0.5*sign(abs(d)+m_hz)*sqrt(1+m_nx/om);
}

std::complex<double> SFpNxpHzWaveFunction::Yk(double* k, size_t band) const
{
    double om=omega(k,band);
    std::complex<double> d=delta(k);
    std::complex<double> pd=d/abs(d);
    if(band==0)
        return 0.5*pd*sqrt(1-m_nx/om);
    if(band==1)
        return -0.5*pd*sqrt(1-m_nx/om);
    if(band==2)
        return 0.5*pd*sqrt(1-m_nx/om);
    else //band==3
        return -0.5*pd*sqrt(1-m_nx/om);
}

std::ostream & operator<<(std::ostream& left, const SFpNxpHzWaveFunction & right)
{
    left<<*((WaveFunction*)(&right))<<std::endl;
    return left;
}
