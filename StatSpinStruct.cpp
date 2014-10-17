#include <cmath>
#include "StatSpinStruct.h"
#include "Stepper.h"
#include "SlaterDeterminant.h"
#include "Jastrow.h"
#include "LatticeState.h"
#include "Lattice.h"

using namespace std;

StatSpinStruct::StatSpinStruct(const Stepper* stepper,
                               FileManager* fm, bool meas_trans)
    :MatrixQuantity(stepper,fm,"StatSpinStruct",5,
                    pow(stepper->GetLatticeState()->GetLattice()->GetNv(),2)), m_meas_trans(meas_trans)
{}

void StatSpinStruct::measure()
{
#ifndef NDEBUG
    cout<<"StatSpinStruct::measure"<<endl;
#endif
    const LatticeState* st=m_stepper->GetLatticeState();
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
    size_t Lx=m_stepper->GetLatticeState()->GetLattice()->GetLx();
    size_t Ly=m_stepper->GetLatticeState()->GetLattice()->GetLy();
    size_t N=Lx*Ly;

    // get spin swap list
    vector<vector<hop_path_t> > hops;
    vector<uint_vec_t> sti,stj;
    uint_vec_t Nifs=st->GetNifs();
    for(size_t vi=0;vi<N;++vi){
        for(size_t vj=vi;vj<N;++vj){
            const Vertex* vxi=st->GetLattice()->GetVertices()[vi];
            const Vertex* vxj=st->GetLattice()->GetVertices()[vj];
            st->GetLatOc(vxi->idx,sti);
            st->GetLatOc(vxj->idx,stj);
            // assuming single occupancy
            size_t fi=max_element(sti.begin(),sti.end(),uint_vec_t_comp)-sti.begin();
            size_t fj=max_element(stj.begin(),stj.end(),uint_vec_t_comp)-stj.begin();
            if(fi!=fj || sti[fi]!=stj[fj]){//+- and -+ components
                hops.push_back(vector<hop_path_t>(st->GetNfl()));
                hops.back()[fi].push_back(hop_t(vxi->idx*Nifs[fi]+sti[fi][0],
                                                vxj->idx*Nifs[fi]+sti[fi][0]));
                hops.back()[fj].push_back(hop_t(vxj->idx*Nifs[fj]+stj[fj][0],
                                                vxi->idx*Nifs[fj]+stj[fj][0]));
            } else if(m_meas_trans && st->GetNfl()==1 && vi!=vj){ //++ and -- components, only if Sztot not conserved.
                if(isup(sti)){
                    hops.push_back(vector<hop_path_t>(st->GetNfl()));
                    hops.back()[0].push_back(hop_t(vxi->idx*2,vxi->idx*2+1));
                    hops.back()[0].push_back(hop_t(vxj->idx*2,vxj->idx*2+1));
                } else {
                    hops.push_back(vector<hop_path_t>(st->GetNfl()));
                    hops.back()[0].push_back(hop_t(vxi->idx*2+1,vxi->idx*2));
                    hops.back()[0].push_back(hop_t(vxj->idx*2+1,vxj->idx*2));
                }
            }
        }
    }

    // Calculate spin swap amplitudes
    vector<BigComplex> swamps(hops.size());
    vector<BigDouble> swjs(hops.size());
    BigComplex amp=m_stepper->GetAmp()->Amp();
    double jas=m_stepper->GetJas()->Jas();
    m_stepper->GetAmp()->VirtUpdate(hops,vector<vector<hop_path_t> >(1,vector<hop_path_t>(st->GetNfl())),swamps);
    m_stepper->GetJas()->VirtUpdate(hops,swjs);

    // Populate sab[i,j] matrices
    complex<double>* szz=&m_vals[0];
    complex<double>* spm=&m_vals[pow(N,2)];
    complex<double>* smp=&m_vals[2*pow(N,2)];
    complex<double>* spp=&m_vals[3*pow(N,2)];
    complex<double>* smm=&m_vals[4*pow(N,2)];
    BigDouble w=m_stepper->weight();
    size_t sw=0;
    for(size_t vi=0; vi < N; ++vi){
        for(size_t vj=vi; vj < N; ++vj){
            const Vertex* vxi=st->GetLattice()->GetVertices()[vi];
            const Vertex* vxj=st->GetLattice()->GetVertices()[vj];
            st->GetLatOc(vxi->idx,sti);
            st->GetLatOc(vxj->idx,stj);
            if(isup(sti)){
                if(isup(stj)){
                    if(vi==vj){
                        spm[vi*N+vj]+=1.0;
                        szz[vi*N+vj]+=0.25;
                    } else {
                        szz[vi*N+vj]+=0.25;
                        szz[vj*N+vi]+=0.25;
                        if(m_meas_trans && st->GetNfl()==1){
                            smm[vi*N+vj]+=complex<double>(conj(amp)*jas*swamps[sw]*swjs[sw]/w);
                            smm[vj*N+vi]+=complex<double>(conj(amp)*jas*swamps[sw]*swjs[sw]/w);
                            ++sw;
                        }
                    }
                } else {
                    szz[vi*N+vj]+=-0.25;
                    szz[vj*N+vi]+=-0.25;
                    smp[vi*N+vj]+=-complex<double>(conj(amp)*jas*swamps[sw]*swjs[sw]/w);
                    spm[vj*N+vi]+=-complex<double>(conj(amp)*jas*swamps[sw]*swjs[sw]/w);
                    ++sw;
                }
            } else {
                if(isup(stj)){
                    szz[vi*N+vj]+=-0.25;
                    szz[vj*N+vi]+=-0.25;
                    spm[vi*N+vj]+=-complex<double>(conj(amp)*jas*swamps[sw]*swjs[sw]/w);
                    smp[vj*N+vi]+=-complex<double>(conj(amp)*jas*swamps[sw]*swjs[sw]/w);
                    ++sw;
                } else {
                    if(vi==vj){
                        smp[vi*N+vj]+=1.0;
                        szz[vi*N+vj]+=0.25;
                    } else {
                        szz[vi*N+vj]+=0.25;
                        szz[vj*N+vi]+=0.25;
                        if(m_meas_trans && st->GetNfl()==1){
                            spp[vi*N+vj]+=complex<double>(conj(amp)*jas*swamps[sw]*swjs[sw]/w);
                            spp[vj*N+vi]+=complex<double>(conj(amp)*jas*swamps[sw]*swjs[sw]/w);
                            ++sw;
                        }
                    }
                }
            }
        }
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

