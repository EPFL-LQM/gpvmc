#include "Magnetization.h"
#include "Stepper.h"
#include "Amplitude.h"
#include "LatticeState.h"
#include "Lattice.h"
#include <array>

using namespace std;

Magnetization::Magnetization(const Stepper* stepper, FileManager* fm)
    :VectorQuantity(stepper,fm,"Magnetization",3)
{
    const LatticeState* st=m_stepper->GetAmp()->GetLatticeState();
    if(!(st->GetNfl()==1 && st->GetNifs()[0]==2))
    {
        string err="Magnetization::Magnatization: only defined for a system "
                   "of spin-1/2 particles with Sztot not conserved.";
#ifdef EXCEPT
        throw(std::logic_error(err.c_str()));
#else
        cerr<<err<<endl;
        abort();
#endif
    }
}

Magnetization::~Magnetization()
{}

void Magnetization::measure()
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
    BigComplex amp=m_stepper->GetAmp()->Amp();
    m_stepper->GetAmp()->VirtUpdate(hops,vector<vector<hop_path_t> >(1,vector<hop_path_t>(1)),swamps);
    vector<BigComplex> S(3,BigComplex(0.0,0.0));
    for(size_t v=0;v<st->GetNsites();++v){
        const Vertex* vxi=st->GetLattice()->GetVertices()[v];
        st->GetLatOc(vxi->idx,sti);
        if(sti[0][0]==0){//up
            S[2]+=0.5;
            S[0]+=0.5*conj(amp)*swamps[v];
        } else {//down
            S[2]-=0.5;
            S[1]+=0.5*conj(amp)*swamps[v];
        }
    }
    BigComplex mI(0.0,-1.0);
    BigDouble w=m_stepper->weight();
    m_vals[0]+=complex<double>((S[0]+S[1])/w)/double(st->GetNsites());
    m_vals[1]+=complex<double>(mI*(S[0]-S[1])/w)/double(st->GetNsites());
    m_vals[2]+=complex<double>(S[2])/double(st->GetNsites());
}
