#include <cmath>
#include "StatSpinStruct.h"
#include "Stepper.h"
#include "SlaterDeterminant.h"
#include "Jastrow.h"
#include "LatticeState.h"
#include "Lattice.h"

using namespace std;

StatSpinStruct::StatSpinStruct(const Stepper* stepper,
                               FileManager* fm)
    :MatrixQuantity(stepper,fm,"StatSpinStruct",3,
                      stepper->GetAmp()->GetLatticeState()->GetLattice()->GetLx()*
                      stepper->GetAmp()->GetLatticeState()->GetLattice()->GetLy())
{
    size_t Lx=stepper->GetAmp()->GetLatticeState()->GetLattice()->GetLx();
    size_t Ly=stepper->GetAmp()->GetLatticeState()->GetLattice()->GetLy();
    m_qs.reserve(Lx*Ly);
    for(size_t qx=0; qx<Lx;++qx){
        for(size_t qy=0;qy<Ly;++qy){
            m_qs.push_back(qx);
            m_qs.push_back(qy);
        }
    }
    complex<double> I(0,1);
    m_ph=vector<complex<double> >(Lx*Ly*m_qs.size()/2);
    for(size_t q=0;q<m_qs.size()/2;++q){
        for(size_t x=0;x<Lx;++x){
            for(size_t y=0;y<Ly;++y){
                double phase=2*M_PI*(double(m_qs[2*q]*x)/Lx+double(m_qs[2*q+1]*y)/Ly);
                m_ph[(q*Lx+x)*Ly+y]=1.0/sqrt(Lx*Ly)*(cos(phase)+I*sin(phase));
            }
        }
    }
}

void StatSpinStruct::measure()
{
    const LatticeState* st=m_stepper->GetAmp()->GetLatticeState();
    if(!(st->GetNfl()==1 && st->GetNifs()[0]==2) &&
                    !(st->GetNfl()==2 && st->GetNifs()[0]==1 && st->GetNifs()[1]))
    {
        string err="StatSpinStruct::measure: Only defined for a system of spin-1/2 particles";
#ifdef EXCEPT
        throw(std::logic_error(err.c_str()));
#else
        cerr<<err<<endl;
        abort();
#endif
    }
    Quantity::measure();
    size_t Lx=m_stepper->GetAmp()->GetLatticeState()->GetLattice()->GetLx();
    size_t Ly=m_stepper->GetAmp()->GetLatticeState()->GetLattice()->GetLy();
    // get spin swap list
    vector<vector<hop_path_t> > hops;
    vector<uint_vec_t> sti,stj;
    uint_vec_t Nifs=st->GetNifs();
    for(size_t vi=0;vi<st->GetNsites();++vi){
        for(size_t vj=vi;vj<st->GetNsites();++vj){
            const Vertex* vxi=st->GetLattice()->GetVertices()[vi];
            const Vertex* vxj=st->GetLattice()->GetVertices()[vj];
            st->GetLatOc(vxi->idx,sti);
            st->GetLatOc(vxj->idx,stj);
            // assuming single occupancy
            size_t fi=max_element(sti.begin(),sti.end(),uint_vec_t_comp)-sti.begin();
            size_t fj=max_element(stj.begin(),stj.end(),uint_vec_t_comp)-stj.begin();
            if(fi!=fj || sti[fi]!=stj[fj]){
                hops.push_back(vector<hop_path_t>(st->GetNfl()));
                hops.back()[fi].push_back(hop_t(vxi->idx*Nifs[fi]+sti[fi][0],
                                                vxj->idx*Nifs[fi]+sti[fi][0]));
                hops.back()[fj].push_back(hop_t(vxj->idx*Nifs[fj]+stj[fj][0],
                                                vxi->idx*Nifs[fj]+stj[fj][0]));
            }
        }
    }
    // Calculate spin swap amplitudes
    vector<BigComplex> swamps(hops.size());
    vector<double> swjs(hops.size());
    BigComplex amp=m_stepper->GetAmp()->Amp();
    double jas=m_stepper->GetJas()->Jas();
    m_stepper->GetAmp()->VirtUpdate(hops,vector<vector<hop_path_t> >(1,vector<hop_path_t>(st->GetNfl())),swamps);
    m_stepper->GetJas()->VirtUpdate(hops,swjs);
    // Calculate Static spin structure factor
    vector<complex<double> > sqlong(m_qs.size()/2,0);
    vector<BigComplex> sqtranspm(m_qs.size()/2,0);
    vector<BigComplex> sqtransmp(m_qs.size()/2,0);
    for(size_t q=0;q<m_qs.size()/2;++q){
        size_t sw=0;
        for(size_t vi=0;vi<st->GetNsites();++vi){
            for(size_t vj=vi;vj<st->GetNsites();++vj){
                const Vertex* vxi=st->GetLattice()->GetVertices()[vi];
                const Vertex* vxj=st->GetLattice()->GetVertices()[vj];
                st->GetLatOc(vxi->idx,sti);
                st->GetLatOc(vxj->idx,stj);
                if(isup(sti)){
                    if(isup(stj)){
                        if(vi==vj){
                            sqtranspm[q]+=norm(amp)*pow(jas,2)/double(Lx*Ly);
                            sqlong[q]+=0.25/double(Lx*Ly);
                        } else {
                            sqlong[q]+=conj(m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]])*
                                            0.25*
                                            m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]];
                            sqlong[q]+=conj(m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]])*
                                            0.25*
                                            m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]];
                        }
                    } else {
                        sqlong[q]+=conj(m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]])*
                                        (-0.25)*
                                        m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]];
                        sqlong[q]+=conj(m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]])*
                                        (-0.25)*
                                        m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]];
                        sqtransmp[q]-=conj(m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]])*
                                           conj(amp)*jas*swamps[sw]*swjs[sw]*
                                           m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]];
                        sqtranspm[q]-=conj(m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]])*
                                           amp*jas*conj(swamps[sw])*swjs[sw]*
                                           m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]];
                        ++sw;
                    }
                } else {
                    if(isup(stj)){
                        sqlong[q]+=conj(m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]])*
                                        (-0.25)*
                                        m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]];
                        sqlong[q]+=conj(m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]])*
                                        (-0.25)*
                                        m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]];
                        sqtranspm[q]-=conj(m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]])*
                                           conj(amp)*jas*swamps[sw]*swjs[sw]*
                                           m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]];
                        sqtransmp[q]-=conj(m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]])*
                                           amp*jas*conj(swamps[sw])*swjs[sw]*
                                           m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]];
                        ++sw;
                    } else {
                        if(vi==vj){
                            sqtransmp[q]+=norm(amp)*pow(jas,2)/double(Lx*Ly);
                            sqlong[q]+=0.25/double(Lx*Ly);
                        } else {
                            sqlong[q]+=conj(m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]])*
                                            0.25*
                                            m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]];
                            sqlong[q]+=conj(m_ph[(q*Lx+vxi->uc[0])*Ly+vxi->uc[1]])*
                                            0.25*
                                            m_ph[(q*Lx+vxj->uc[0])*Ly+vxj->uc[1]];
                        }
                    }
                }
            }
        }
    }
    BigDouble w=m_stepper->weight();
    for(size_t q=0;q<m_qs.size()/2;++q){
        Val(0,q)+=sqlong[q];
        Val(1,q)+=complex<double>(sqtransmp[q]/w);
        Val(2,q)+=complex<double>(sqtranspm[q]/w);
    }
}

bool StatSpinStruct::isup(const vector<uint_vec_t>& in)
{
    if(in.size()==1)
        return in[0].size()==1 && in[0][0]==0;
    else
        return in[0].size()==1 && in[1].size()==0;
}

bool StatSpinStruct::isdo(const vector<uint_vec_t>& in)
{
    if(in.size()==1)
        return in[0].size()==1 && in[0][0]==1;
    else
        return in[0].size()==0 && in[1].size()==1;
}

