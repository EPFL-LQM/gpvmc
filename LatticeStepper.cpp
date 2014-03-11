#include "LatticeStepper.h"
#include "Amplitude.h"
#include "LatticeState.h"
#include "Lattice.h"
#include "WaveFunction.h"
#include "linalg.h"
#include "RanGen.h"
#include "Timer.h"
#include <numeric>

using namespace std;

LatticeStepper::LatticeStepper(Amplitude* amp)
    :Stepper(amp),
     m_Nfl(amp->GetLatticeState()->GetNfl()),
     m_Nifs(amp->GetLatticeState()->GetNifs()),
     m_prev(amp->GetLatticeState()->GetNfl()),
     m_khop(-1), m_weight(-1)
{}


void LatticeStepper::Reset()
{
    m_prev=vector<hop_path_t>(m_Nfl);
    m_khop=-1;
    m_prev_weight=-1;
    m_weight=-1;
}

BigDouble LatticeStepper::trystep()
{
#ifdef PROFILE
    Timer::tic("LatticeStepper::trystep");
#endif
    const WaveFunction* wav=m_amp->GetWaveFunction();
    const LatticeState* st=m_amp->GetLatticeState();
    const Lattice* lat=st->GetLattice();
    const Vertex* v1,*v2;
    // it is assumed all vertices hold exactly one fermion
    // pick up a random vertex
    v1=lat->GetVertices()[size_t(lat->GetNv()*RanGen::uniform())];
    // pick up a neighbour at random
    v2=v1->edges[size_t(v1->edges.size()*RanGen::uniform())]->GetOther(v1);
    vector<uint_vec_t> s1,s2;
    st->GetLatOc(v1->idx,s1);
    st->GetLatOc(v2->idx,s2);
    // get s1 occupied flavour
    size_t f1=max_element(s1.begin(),s1.end(),
                          [](const uint_vec_t& a, const uint_vec_t& b)
                          {return a.size()<b.size();})-s1.begin();
    // get s2 occupied flavour
    size_t f2=max_element(s2.begin(),s2.end(),
                          [](const uint_vec_t& a, const uint_vec_t& b)
                          {return a.size()<b.size();})-s2.begin();
    // from here it is assumed there is only one fermion per site.
    size_t c= m_Nifs[f1]*RanGen::uniform();
    if(c==s1[f1][0]){
        // exchange neighbouring fermions
        m_prev=vector<hop_path_t>(m_Nfl);
        if(f1==f2 && s1[f1][0]==s2[f2][0]){
            // do nothing if neighbouring fermions are in same state
            m_khop=-1;
            m_prev_weight=m_weight;
            return m_weight;
        } else {
            m_prev[f1].push_back(hop_t(v1->idx*m_Nifs[f1]+s1[f1][0],
                                       v2->idx*m_Nifs[f1]+s1[f1][0]));
            m_prev[f2].push_back(hop_t(v2->idx*m_Nifs[f2]+s2[f2][0],
                                       v1->idx*m_Nifs[f2]+s2[f2][0]));
        }
    } else {
        // flip fermion on v1 to any of the
        // other internal degrees of freedom
        m_prev=vector<hop_path_t>(m_Nfl);
        m_prev[f1].push_back(hop_t(v1->idx*m_Nifs[f1]+s1[f1][0],v1->idx*m_Nifs[f1]+c));
    }
    const vector<vector<hop_path_t> >& khops=wav->GetHops();
    vector<BigComplex> qs(khops.size());
#ifdef PROFILE
    Timer::tic("LatticeStepper::trystep/VirtUpdate");
#endif
    m_amp->VirtUpdate(vector<vector<hop_path_t> >(1,m_prev),khops,qs);
#ifdef PROFILE
    Timer::toc("LatticeStepper::trystep/VirtUpdate");
#endif
    BigDouble out(0);
    for_each(qs.begin(),qs.end(),[&](const BigComplex& aq){out+=norm(aq);}); 
    m_khop=max_element(qs.begin(),qs.end(),
                       [](const BigComplex& a,const BigComplex& b)
                       {return norm(a)<norm(b);}) - qs.begin();
    m_prev_weight=out;
#ifdef PROFILE
    Timer::toc("LatticeStepper::trystep");
#endif
    return out;
}

void LatticeStepper::step()
{
    LatticeState* st=m_amp->GetLatticeState();
    WaveFunction *wav=m_amp->GetWaveFunction();
    st->Hop(m_prev);
    if(!m_khop<0){
        const vector<hop_path_t>& khop_p=wav->GetHop(m_khop);
        wav->Hop(m_khop);
        m_amp->Update(m_prev,khop_p);
    } else {
        m_amp->Update(m_prev,vector<hop_path_t>(m_Nfl));
    }
    if(m_amp->Amp()==0){
        cerr<<"LatticeStepper::step(): Stepped to a zero overlap!"<<endl;
        abort();
    }
    BigDouble ow=m_weight;
#ifndef DEBUG
    if(m_prev_weight>0.0)
        m_weight=m_prev_weight;
    else {
        m_weight=-1;
        weight();
    }
#else
    m_weight=-1;
    weight();
    if(m_prev_weight>=0.0 && abs(double(m_prev_weight/m_weight)-1)>1e-6){
        cerr<<"LatticeStepper::step: error: abs(double(m_prev_weight/m_weight))="<<abs(double(m_prev_weight/m_weight))-1<<" m_prev_weight="<<m_prev_weight<<" and m_weight="<<m_weight<<endl;
        abort();
    }
#endif
    if(abs(double(m_weight/ow))<1e-6){
        cerr<<"LatticeStepper::step: warning, abs(double(m_weight/ow)) is VERY small! : "<<abs(double(m_weight/ow))<<endl;
    }
    m_prev=vector<hop_path_t>(m_Nfl);
    m_khop=-1;
}

BigDouble LatticeStepper::weight()
{
    const WaveFunction* wav=m_amp->GetWaveFunction();
    if(m_weight<0.0){
        BigDouble out(0);
        const vector<vector<hop_path_t> >& khop=wav->GetHops();
        vector<vector<hop_path_t> > rhop(1,vector<hop_path_t>(2));
        vector<BigComplex> qs(khop.size(),0);
        m_amp->VirtUpdate(rhop,khop,qs);
        for(size_t k=0;k<khop.size();++k){
            out+=norm(qs[k]);
        }
        m_weight=out;
        if(m_weight==0.0){
            std::cout<<"Warning: LatticeStepper::weight: "
                  "state has no overlap."<<std::endl;
        }
        return m_weight;
    } else {
        return m_weight;
    }
}

double LatticeStepper::transprob()
{
    return 1.0;
}
