#include "SFpNpHxWaveFunction.h"
#include <cstdlib>
#include <stdexcept>
#include <iostream>

using namespace std;

SFpNpHxWaveFunction::SFpNpHxWaveFunction(FileManager* fm,
                                         size_t Lx, size_t Ly,
                                         size_t Nby,
                                         double phi, double neel, double hx,
                                         const vector<double>& bc_phase)
    :WaveFunction(fm),
     m_Lx(Lx), m_Ly(Ly),
     m_qn2fock(new size_t[Lx*Ly*4]),
     m_fock2qn(new size_t[Lx*Ly*2*3]),
     m_phi(phi), m_neel(neel), m_hx(hx),
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

SFpNpHxWaveFunction::~SFpNpHxWaveFunction()
{
    delete [] m_fock2qn;
    delete [] m_qn2fock;
}

std::complex<double> SFpNpHxWaveFunction::matrix_element(size_t fk,size_t fr, size_t fl)
{
    std::complex<double> I(0,1);
    double k[2], r[2];
    k[0]=(m_fock2qn[3*fk]+m_bc_phase[0]/(2.0))*2*M_PI/double(m_Lx);
    k[1]=(m_fock2qn[3*fk+1]+m_bc_phase[1]/(2.0))*2*M_PI/double(m_Ly);
    int ifr=fr;
    size_t up=ifr%2;
    r[1]=((ifr-up)/2)%m_Ly;
    r[0]=((ifr-up-r[1]*2)/(2*m_Ly));
    up=size_t(1-int(up));
    size_t band=m_fock2qn[3*fk+2];
    std::complex<double> out;
    out= sqrt(2.0)*exp(I*(k[0]*r[0]+k[1]*r[1]))*
                (double((1-int(r[0]+r[1])%2))*conj(Uk(k,up,band))
                +double(int(r[0]+r[1])%2)*conj(Vk(k,up,band)));
#ifdef CRAY
    if(isnan(real(out)) || isnan(imag(out)))
#else
    if(std::isnan(real(out)) || std::isnan(imag(out)))
#endif
    {
#ifdef EXCEPT
        throw std::logic_error("SFpNpHxWaveFunction::matrix_element::NaN encoutered");
#else
        std::cerr<<"SFpNpHxWaveFunction::matrix_element::NaN encoutered"<<std::endl;
        abort();
#endif
    }
    return out;
}

void SFpNpHxWaveFunction::quantum_numbers(const size_t& f, const size_t& fl, map<string,size_t>& qn)
{
    qn.clear();
    qn["kx"]=m_fock2qn[3*f];
    qn["ky"]=m_fock2qn[3*f+1];
    qn["band"]=m_fock2qn[3*f+2];
}

bool SFpNpHxWaveFunction::inmbz(size_t kx, size_t ky) const
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

void SFpNpHxWaveFunction::mbzmod(size_t* k) const
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

int SFpNpHxWaveFunction::sign(double x) const
{
    if(x>0)
        return 1;
    else if(x<0)
        return -1;
    else
        return 0;
}

std::complex<double> SFpNpHxWaveFunction::delta(double* k) const
{
    return 0.5*(cos(k[0])*std::polar(1.0,m_phi)+cos(k[1])*std::polar(1.0,-m_phi));
}

double SFpNpHxWaveFunction::omega(double* k, size_t band) const
{
    if(band==0){
        return -sqrt(m_neel*m_neel + pow(m_hx-abs(delta(k)),2));
    } else if(band==1){
        return -sqrt(m_neel*m_neel + pow(m_hx+abs(delta(k)),2));
    } else if(band==2){
        return sqrt(m_neel*m_neel + pow(m_hx-abs(delta(k)),2));
    } else {//band==3
        return sqrt(m_neel*m_neel + pow(m_hx+abs(delta(k)),2));
    }
}

std::complex<double> SFpNpHxWaveFunction::Uk(double* k, bool up, size_t band) const
{
    double om=omega(k,band);
    std::complex<double> d=delta(k);
    if(band==0){
        if(up){
            return 0.5*sqrt(1.0-m_neel/om);
        } else {
            return 0.5*sign(abs(d)-m_hx)*sqrt(1.0+m_neel/om);
        }
    } else if(band==1){
        if(up){
            return 0.5*sqrt(1.0-m_neel/om);
        } else {
            return -0.5*sign(abs(d)+m_hx)*sqrt(1.0+m_neel/om);
        }
    } else if(band==2){
        if(up){
            return 0.5*sqrt(1.0-m_neel/om);
        } else {
            return -0.5*sign(abs(d)-m_hx)*sqrt(1.0+m_neel/om);
        }
    } else {//band==3
        if(up){
            return 0.5*sqrt(1.0-m_neel/om);
        } else {
            return 0.5*sign(abs(d)+m_hx)*sqrt(1.0+m_neel/om);
        }
    }
}

std::complex<double> SFpNpHxWaveFunction::Vk(double* k, bool up, size_t band) const
{
    double om=omega(k,band);
    std::complex<double> d=delta(k);
    complex<double> pd=conj(d)/abs(d);
    if(band==0){
        if(up){
            return 0.5*sign(abs(d)-m_hx)*pd*sqrt(1.0+m_neel/om);
        } else {
            return 0.5*pd*sqrt(1.0-m_neel/om);
        }
    } else if(band==1){
        if(up){
            return 0.5*sign(abs(d)+m_hx)*pd*sqrt(1.0+m_neel/om);
        } else {
            return -0.5*pd*sqrt(1.0-m_neel/om);
        }
    } else if(band==2){
        if(up){
            return -0.5*sign(abs(d)-m_hx)*pd*sqrt(1.0+m_neel/om);
        } else {
            return 0.5*pd*sqrt(1.0-m_neel/om);
        }
    } else {//band==3
        if(up){
            return -0.5*sign(abs(d)+m_hx)*pd*sqrt(1.0+m_neel/om);
        } else {
            return -0.5*pd*sqrt(1.0-m_neel/om);
        }
    }
}

std::ostream & operator<<(std::ostream& left, const SFpNpHxWaveFunction & right)
{
    left<<*((WaveFunction*)(&right))<<std::endl;
    return left;
}
