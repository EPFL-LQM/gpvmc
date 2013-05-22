#include "StagFluxWaveFunction.h"
#include <cstdlib>
#include <stdexcept>
#include <iostream>

StagFluxWaveFunction::StagFluxWaveFunction(size_t Lx, size_t Ly,
                                           size_t Nfsup, size_t Nfsdo,
                                           double phi, double neel, double *bc_phase)
    :WaveFunction(Lx,Ly,Nfsup,Nfsdo),m_qn2fock(new size_t[Lx*Ly*2]),
     m_fock2qn(new size_t[Lx*Ly*3]), m_phi(phi), m_neel(neel)
{
    m_bc_phase[0]=bc_phase[0];
    m_bc_phase[1]=bc_phase[1];
    size_t idx=0;
    for(size_t i=0;i<2*Lx*Ly;++i) m_qn2fock[i]=Lx*Ly;
    for(size_t kx=0;kx<Lx;++kx){
        for(size_t ky=0;ky<Ly;++ky){
            if(inmbz(kx,ky)){
                m_fock2qn[idx*3]=kx;
                m_fock2qn[idx*3+1]=ky;
                m_fock2qn[idx*3+2]=0;
                m_qn2fock[(kx*Ly+ky)*2]=idx;
                idx+=1;
                m_fock2qn[idx*3]=kx;
                m_fock2qn[idx*3+1]=ky;
                m_fock2qn[idx*3+2]=1;
                m_qn2fock[(kx*Ly+ky)*2+1]=idx;
                idx+=1;
            }
        }
    }
    init_matrices();
}

StagFluxWaveFunction::~StagFluxWaveFunction()
{
    delete [] m_fock2qn;
    delete [] m_qn2fock;
}

std::complex<double> StagFluxWaveFunction::matrix_element(size_t fk,
                                                     size_t fr, bool up)
{
    std::complex<double> I(0,1);
    double k[2], r[2];
    k[0]=(m_fock2qn[3*fk]+m_bc_phase[0]/(2.0))*2*M_PI/double(m_Lx);
    k[1]=(m_fock2qn[3*fk+1]+m_bc_phase[1]/(2.0))*2*M_PI/double(m_Ly);
    r[1]=fr/m_Lx;
    r[0]=fr-r[1]*m_Lx;
    size_t band=m_fock2qn[3*fk+2];
    std::complex<double> out;
    out= sqrt(2.0)*exp(I*(k[0]*r[0]+k[1]*r[1]))*
                (double((1-int(r[0]+r[1])%2))*Uk(k,up,band)
                +double(int(r[0]+r[1])%2)*Vk(k,up,band));
#ifdef CRAY
    if(isnan(real(out)) || isnan(imag(out)))
#else
    if(std::isnan(real(out)) || std::isnan(imag(out)))
#endif
#ifdef EXCEPT
        throw std::logic_error("StagFluxWaveFunction::matrix_element::NaN encoutered");
#else
        std::cerr<<"StagFluxWaveFunction::matrix_element::NaN encoutered"<<std::endl;
        abort();
#endif
    return out;
}

bool StagFluxWaveFunction::inmbz(size_t kx, size_t ky) const
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

void StagFluxWaveFunction::mbzmod(size_t* k) const
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

std::complex<double> StagFluxWaveFunction::delta(double* k) const
{
    return 0.5*(cos(k[0])*std::polar(1.0,m_phi)+cos(k[1])*std::polar(1.0,-m_phi));
}

double StagFluxWaveFunction::omega(double* k) const
{
    return sqrt(m_neel*m_neel + norm(delta(k)));
}

std::complex<double> StagFluxWaveFunction::Uk(double* k, bool up, size_t band) const
{
    double om=omega(k);
    std::complex<double> d=delta(k);
    int sig=2*int(up)-1;
    if(band==0){
        return sqrt(0.5*std::complex<double>(1+sig*m_neel/om));
    } else {
        return -conj(d)/abs(d)*sqrt(0.5*std::complex<double>(1-sig*m_neel/om));
    }
}

std::complex<double> StagFluxWaveFunction::Vk(double* k, bool up, size_t band) const
{
    double om=omega(k);
    std::complex<double> d=delta(k);
    int sig=2*int(up)-1;
    if(band==0){
        return d/abs(d)*sqrt(0.5*std::complex<double>(1-sig*m_neel/om));
    } else {
        return sqrt(0.5*std::complex<double>(1+sig*m_neel/om));
    }
}

std::ostream & operator<<(std::ostream& left, const StagFluxWaveFunction & right)
{
    left<<*((WaveFunction*)(&right))<<std::endl;
    left<<"spins up quantum number and fock space index:"<<std::endl;
    for(size_t f=0;f<right.m_Nfsup;++f)
        left<<std::setw(2)<<f<<" ";
    left<<std::endl<<"------------------------------------------------"<<std::endl;
    for(size_t f=0;f<right.m_Nfsup;++f)
        left<<std::setw(2)<<right.m_fsup[f]<<" ";
    left<<std::endl;
    for(size_t f=0;f<right.m_Nfsup;++f)
        left<<std::setw(2)<<right.m_fock2qn[3*right.m_fsup[f]]<<" ";
    left<<std::endl;
    for(size_t f=0;f<right.m_Nfsup;++f)
        left<<std::setw(2)<<right.m_fock2qn[3*right.m_fsup[f]+1]<<" ";
    left<<std::endl;
    for(size_t f=0;f<right.m_Nfsup;++f)
        left<<std::setw(2)<<right.m_fock2qn[3*right.m_fsup[f]+2]<<" ";
    left<<std::endl;
    left<<"spins do quantum number and fock space index:"<<std::endl;
    for(size_t f=0;f<right.m_Nfsdo;++f)
        left<<std::setw(2)<<f<<" ";
    left<<std::endl<<"------------------------------------------------"<<std::endl;
    for(size_t f=0;f<right.m_Nfsdo;++f)
        left<<std::setw(2)<<right.m_fsdo[f]<<" ";
    left<<std::endl;
    for(size_t f=0;f<right.m_Nfsdo;++f)
        left<<std::setw(2)<<right.m_fock2qn[3*right.m_fsdo[f]]<<" ";
    left<<std::endl;
    for(size_t f=0;f<right.m_Nfsdo;++f)
        left<<std::setw(2)<<right.m_fock2qn[3*right.m_fsdo[f]+1]<<" ";
    left<<std::endl;
    for(size_t f=0;f<right.m_Nfsdo;++f)
        left<<std::setw(2)<<right.m_fock2qn[3*right.m_fsdo[f]+2]<<" ";
    left<<std::endl;
    return left;
}
