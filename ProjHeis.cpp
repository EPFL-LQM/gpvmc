#include "ProjHeis.h"
#ifdef USEPARA
#include <omp.h>
#endif

#include "Stepper.h"
#include "Jastrow.h"
#include "linalg.h"
#include <vector>
#include <iostream>

using namespace std;

ProjHeis::ProjHeis(const Stepper* stepper,
                   FileManager* fm,
                   double jr, double hz)
    :MatrixQuantity(stepper,fm,"ProjHeis",
            2*stepper->GetAmp()->GetWaveFunction()->GetNExc(),
            stepper->GetAmp()->GetWaveFunction()->GetNExc()),
        m_jr(jr), m_hz(hz)
{}

void ProjHeis::measure()
{
#ifdef PROFILE
    Timer::tic("ProjHeis::measure");
#endif
    Quantity::measure();
    const Amplitude* amp=m_stepper->GetAmp();
    const SpinState* st=amp->GetSpinState();
    const WaveFunction* wav=amp->GetWaveFunction();
    size_t Nexc=wav->GetNExc();
    BigDouble weight=m_stepper->weight();
    BigComplex *amps=new BigComplex[2*Nexc];
    BigComplex *heisamps=&amps[Nexc];
    int taus[8]={1,0,0,1,1,1,1,-1};//,-1,0,0,-1};
    double Js[4]={1,1,m_jr,m_jr};
    size_t nm=2;
    if(m_jr!=0) nm=4;
    vector<hop_path_t> rhopup, rhopdo,khopup,khopdo;
    wav->GetHops(khopup,khopdo);
    std::vector<double> vJs;
    for(size_t k=0;k<Nexc;++k){
        heisamps[k]=0;
        amps[k]=0;
    }
    // create real space swap hops.
    for(size_t ix=0;ix<st->GetL();++ix){
        for(size_t iy=0;iy<st->GetL();++iy){
            for(size_t t=0;t<nm;++t){
                size_t jx=linalg::mod(int(ix)+taus[2*t],
                                      st->GetL());
                size_t jy=linalg::mod(int(iy)+taus[2*t+1],
                                      st->GetL());
                if(st->GetLatOc(ix,iy)!=st->GetLatOc(jx,jy)){
                    if(st->GetLatOc(ix,iy)==UP){
                        rhopup.push_back(
                                hop_path_t(1,
                                    pair<size_t,size_t>(
                                        st->GetLatUpId(ix,iy),
                                        jy*st->GetL()+jx)
                                    )
                                );
                        rhopdo.push_back(
                                hop_path_t(1,
                                    pair<size_t,size_t>(
                                        st->GetLatDoId(jx,jy),
                                        iy*st->GetL()+ix)
                                    )
                                );
                    } else {
                        rhopup.push_back(
                                hop_path_t(1,
                                    pair<size_t,size_t>(
                                        st->GetLatUpId(jx,jy),
                                        iy*st->GetL()+ix)
                                    )
                                );
                        rhopdo.push_back(
                                hop_path_t(1,
                                    pair<size_t,size_t>(
                                        st->GetLatDoId(ix,iy),
                                        jy*st->GetL()+jx)
                                    )
                                );
                    }
                    vJs.push_back(Js[t]);
                }
            }
        }
    }
    // add the no swap empty hop path
    rhopup.push_back(hop_path_t());
    rhopdo.push_back(hop_path_t());
    vector<BigComplex> qs(rhopup.size()*Nexc);
#ifdef PROFILE
    Timer::tic("ProjHeis::measure/VirtUpdate");
#endif
    amp->VirtUpdate(rhopup,rhopdo,khopup,khopdo,&qs[0]);
    // qs has rhopup.size() lines and khopup.size() columns
#ifdef PROFILE
    Timer::toc("ProjHeis::measure/VirtUpdate");
#endif
    // now fill the amplitudes:
    // amps[k]=<a|k>
    // heisamps[k]=sum_b sum_ij H_ij <b|k>
    // sum over b is trunkated to sum over
    // the |b> defined by hops (hopup,hopdo).
    size_t Nsw=rhopup.size();
    for(size_t k=0;k<Nexc;++k){
        amps[k]=qs[k*Nsw+Nsw-1];
        if(st->GetJas()) amps[k]*=st->GetJas()->Jas();
    }
    size_t swc=0;
    for(size_t ix=0;ix<st->GetL();++ix){
        for(size_t iy=0;iy<st->GetL();++iy){
	    if(st->GetLatOc(ix,iy)==UP){
                for(size_t k=0;k<Nexc;++k)
                    heisamps[k]-=m_hz*amps[k];
            } else if(st->GetLatOc(ix,iy)==DOWN){
                for(size_t k=0;k<Nexc;++k)
                    heisamps[k]+=m_hz*amps[k];
            }
            for(size_t t=0;t<nm;++t){
                size_t jx=linalg::mod(int(ix)+taus[2*t],
                                      st->GetL());
                size_t jy=linalg::mod(int(iy)+taus[2*t+1],
                                      st->GetL());
                if(st->GetLatOc(ix,iy)==st->GetLatOc(jx,jy)){
                    for(size_t k=0;k<Nexc;++k)
                        heisamps[k]+=Js[t]*0.25*amps[k];
                } else {
                    if(st->GetJas()){
                        for(size_t k=0;k<Nexc;++k){
                            heisamps[k]-=0.25*Js[t]*amps[k]+
                                         0.5*Js[t]*qs[k*Nsw+swc]*
                                         st->GetJas()->virtualhop(
                                                 rhopup[swc],
                                                 rhopdo[swc]);
                        }
                    } else {
                        for(size_t k=0;k<Nexc;++k){
                            heisamps[k]-=0.25*Js[t]*amps[k]+
                                         0.5*Js[t]*qs[k*Nsw+swc];
                        }
                    }
                    ++swc;
                }
            }
        }
    }
    // now fill the H and O matrices
    if(Nexc==1){
        Val(0,0)+=std::complex<double>(heisamps[0]*conj(amps[0])/weight/st->GetL()/st->GetL());
        Val(0,1)+=1;
    } else {
        for(size_t k1=0;k1<Nexc;++k1){
            for(size_t k2=0;k2<Nexc;++k2){
                Val(k1,k2)+=std::complex<double>(
                        conj(amps[k1])*heisamps[k2]/weight);
                Val(Nexc+k1,k2)+=std::complex<double>(
                        conj(amps[k1])*amps[k2]/weight);
            }
        }
    }
    delete [] amps;
#ifdef PROFILE
    Timer::toc("ProjHeis::measure");
#endif
}

std::string ProjHeis::str() const
{
    size_t Nks=m_stepper->GetAmp()->GetWaveFunction()->GetNExc();
    std::ostringstream sout;
    sout<<std::scientific<<std::setprecision(10);
    sout<<"Hkk="<<std::endl
        <<linalg::PrintMat(&Val(0,0),Nks,Nks,0,2,false)
        <<std::endl
        <<"Okk="<<std::endl
        <<linalg::PrintMat(&Val(Nks,0),Nks,Nks,0,2,false)
        <<std::endl;
    return sout.str();
}
