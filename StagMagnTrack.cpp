#include "StagMagnTrack.h"
#include "Stepper.h"
#include "Amplitude.h"
#include "LatticeState.h"
#include "FileManager.h"
#include "Lattice.h"

StagMagnTrack::StagMagnTrack(const Stepper* stepper, FileManager* fm)
    :Quantity(stepper,fm,"StagMagnTrack")
{}
StagMagnTrack::~StagMagnTrack() {}
void StagMagnTrack::init()
{
    Quantity::init();
    m_vals.clear();
}
void StagMagnTrack::save() const
{
    m_fm->FileStream(m_basename).Write(m_vals.size()*2,1,(double*)m_vals.data());
}
void StagMagnTrack::measure()
{
    const LatticeState* st=m_stepper->GetAmp()->GetLatticeState();
    if((st->GetNfl()==1 && st->GetNifs()[0]==2) ||
            (st->GetNfl()==2 && st->GetNifs()[0]==1 && st->GetNifs()[1]))
    {
        m_vals.push_back(0);
        const Lattice* lat=st->GetLattice();
        for(size_t v=0;v<st->GetNsites();++v){
            complex<double> ph=std::exp(
                    complex<double>(0,1)*M_PI*(
                        lat->GetVertices()[v]->uc[0]
                        +lat->GetVertices()[v]->pos[0]
                        +lat->GetVertices()[v]->uc[0]
                        +lat->GetVertices()[v]->pos[0]));
            vector<uint_vec_t> vst;
            st->GetLatOc(v,vst);
            int up=0;
            if(st->GetNfl()==2){
                if(vst[0].size())
                    up=1;
            } else {
                if(vst[0][0]==0)
                    up=1;
            }
            m_vals.back()+=double(2*up-1)*ph;
        }
        m_vals.back()/=st->GetNsites();
    }
}


