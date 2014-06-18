#include "StagMagn.h"
#include "Stepper.h"
#include "SlaterDeterminant.h"
#include "Jastrow.h"
#include "LatticeState.h"
#include "Lattice.h"

StagMagn::StagMagn(const Stepper* stepper, FileManager* fm)
    :VectorQuantity(stepper,fm,"StagMagn",3)
{}
StagMagn::~StagMagn()
{}
void StagMagn::measure()
{
    const LatticeState* st=m_stepper->GetAmp()->GetLatticeState();    
    Quantity::measure();
    vector<vector<hop_path_t> > hops;
    vector<uint_vec_t> sti;
    for(size_t v=0;v<st->GetNsites();++v){
        const Vertex* vxi=st->GetLattice()->GetVertices()[v];
        st->GetLatOc(vxi->idx,sti);
        size_t fi=sti[0][0];
        if(fi==0){//up
            hops.push_back(vector<hop_path_t>(1));
            hops.back()[0].push_back(hop_t(vxi->idx*2,vxi->idx*2+1));
        } else {
            hops.push_back(vector<hop_path_t>(1));
            hops.back()[0].push_back(hop_t(vxi->idx*2+1,vxi->idx*2));
        }
    }
    vector<BigComplex> swamps(hops.size(),BigComplex(0.0,0.0));
    vector<double> swjs(hops.size(),0);
    BigComplex amp=m_stepper->GetAmp()->Amp()*m_stepper->GetJas()->Jas();
    m_stepper->GetAmp()->VirtUpdate(hops,vector<vector<hop_path_t> >(1,vector<hop_path_t>(1)),swamps);
    m_stepper->GetJas()->VirtUpdate(hops,swjs);
    vector<BigComplex> S(3,BigComplex(0.0,0.0));
    for(size_t v=0;v<st->GetNsites();++v){
        const Vertex* vxi=st->GetLattice()->GetVertices()[v];
            complex<double> ph=std::exp(
                    complex<double>(0,1)*M_PI*(vxi->uc[0]+vxi->pos[0]+
                                               vxi->uc[1]+vxi->pos[1]));
        st->GetLatOc(vxi->idx,sti);
        if(sti[0][0]==0){//up
            S[2]+=0.5*ph;
            S[0]+=0.5*conj(amp)*swamps[v]*swjs[v]*ph;
        } else {//down
            S[2]-=0.5*ph;
            S[1]+=0.5*conj(amp)*swamps[v]*swjs[v]*ph;
        }
    }
    BigComplex mI(0.0,-1.0);
    BigDouble w=m_stepper->weight();
    m_vals[0]+=complex<double>((S[0]+S[1])/w)/double(st->GetNsites());
    m_vals[1]+=complex<double>(mI*(S[0]-S[1])/w)/double(st->GetNsites());
    m_vals[2]+=complex<double>(S[2])/double(st->GetNsites());
}


