#include "ProjHeis.h"
#ifdef USEPARA
#include <omp.h>
#endif

#include "Stepper.h"
#include "WaveFunction.h"
#include "LatticeState.h"
#include "Amplitude.h"
#include "linalg.h"
#include "Lattice.h"
#include <vector>
#include <iostream>

using namespace std;

ProjHeis::ProjHeis(const Stepper* stepper,
                   FileManager* fm, const Lattice* lat,
                   double Bx)
    :MatrixQuantity(stepper,fm,"ProjHeis",
            2*stepper->GetAmp()->GetWaveFunction()->GetNExc(),
            stepper->GetAmp()->GetWaveFunction()->GetNExc()),
     m_lat(lat), m_Bx(Bx)
{}

void ProjHeis::measure()
{
#ifdef PROFILE
    Timer::tic("ProjHeis::measure");
#endif
    Quantity::measure();
    const Amplitude* amp=m_stepper->GetAmp();
    const LatticeState* st=amp->GetLatticeState();
    const WaveFunction* wav=amp->GetWaveFunction();
    size_t Nexc=wav->GetNExc();
    BigDouble weight=m_stepper->weight();
    vector<BigComplex> amps(Nexc,0),heisamps(Nexc,0);
    vector<vector<hop_path_t> > rhop;
    const vector<vector<hop_path_t> >& khop=wav->GetHops();
    // create real space swap hops.
    uint_vec_t Nifs=st->GetNifs();
    for(size_t e=0;e<m_lat->GetEdges().size();++e){
        vector<uint_vec_t> sti,stj;
        const Vertex* vi=m_lat->GetEdges()[e]->first;
        const Vertex* vj=m_lat->GetEdges()[e]->second;
        st->GetLatOc(vi->idx,sti);
        st->GetLatOc(vj->idx,stj);
        // assuming single occupancy
        size_t fi=max_element(sti.begin(),sti.end(),uint_vec_t_comp)-sti.begin();
        size_t fj=max_element(stj.begin(),stj.end(),uint_vec_t_comp)-stj.begin();
        if(fi!=fj || sti[fi]!=stj[fj]){
            rhop.push_back(vector<hop_path_t>(st->GetNfl()));
            rhop.back()[fi].push_back(hop_t(vi->idx*Nifs[fi]+sti[fi][0],
                                           vj->idx*Nifs[fi]+sti[fi][0]));
            rhop.back()[fj].push_back(hop_t(vj->idx*Nifs[fj]+stj[fj][0],
                                            vi->idx*Nifs[fj]+stj[fj][0]));
        }
    }
    // add the no swap empty hop path
    rhop.push_back(vector<hop_path_t>(st->GetNfl()));
    vector<BigComplex> qs(rhop.size()*Nexc);
#ifdef PROFILE
    Timer::tic("ProjHeis::measure/VirtUpdate");
#endif
    amp->VirtUpdate(rhop,khop,qs);
    // qs has rhop.size() lines and khop.size() columns
#ifdef PROFILE
    Timer::toc("ProjHeis::measure/VirtUpdate");
#endif
    // now fill the amplitudes:
    // amps[k]=<a|k>
    // heisamps[k]=sum_b sum_ij H_ij <b|k>
    // sum over b is trunkated to sum over
    // the |b> defined by hops (hopup,hopdo).
    size_t Nsw=rhop.size();
    for(size_t k=0;k<Nexc;++k){
        amps[k]=qs[k*Nsw+Nsw-1];
    }
    size_t swc=0;
    for(size_t e=0;e<m_lat->GetEdges().size();++e){
        double J=m_lat->GetEdges()[e]->Jprop;
        vector<uint_vec_t> sti,stj;
        const Vertex* vi=m_lat->GetEdges()[e]->first;
        const Vertex* vj=m_lat->GetEdges()[e]->second;
        st->GetLatOc(vi->idx,sti);
        st->GetLatOc(vj->idx,stj);
        // assuming single occupancy
        size_t fi=max_element(sti.begin(),sti.end(),uint_vec_t_comp)-sti.begin();
        size_t fj=max_element(stj.begin(),stj.end(),uint_vec_t_comp)-stj.begin();
        if(fi!=fj || sti[fi]!=stj[fj]){
            for(size_t k=0;k<Nexc;++k)
                heisamps[k]-=J*(0.25*amps[k]+0.5*qs[k*Nsw+swc]);
            ++swc;
        } else {
            for(size_t k=0;k<Nexc;++k)
                heisamps[k]+=m_lat->GetEdges()[e]->Jprop*0.25*amps[k];
        }
    }
    // now fill the H and O matrices
    if(Nexc==1){
        Val(0,0)+=std::complex<double>(heisamps[0]*conj(amps[0])/weight/m_lat->GetLx()/m_lat->GetLy());
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
