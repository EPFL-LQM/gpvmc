#include "StagMagnZ.h"
#include "Stepper.h"
#include "SlaterDeterminant.h"
#include "LatticeState.h"
#include "Lattice.h"

StagMagnZ::StagMagnZ(const Stepper* stepper, FileManager* fm)
    :ScalarQuantity(stepper,fm,"StagMagnZ")
{}
StagMagnZ::~StagMagnZ()
{}
void StagMagnZ::measure()
{
    Quantity::measure();
    const LatticeState* st=m_stepper->GetLatticeState();
    complex<double> msz(0);
    if((st->GetNfl()==1 && st->GetNifs()[0]==2) ||
            (st->GetNfl()==2 && st->GetNifs()[0]==1 && st->GetNifs()[1]==1))
    {
        const Lattice* lat=st->GetLattice();
        for(size_t v=0;v<st->GetNsites();++v){
            complex<double> ph=std::exp(
                    complex<double>(0,1)*M_PI*(
                        lat->GetVertices()[v]->uc[0]
                        +lat->GetVertices()[v]->pos[0]
                        +lat->GetVertices()[v]->uc[1]
                        +lat->GetVertices()[v]->pos[1]));
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
            msz+=ph*double(2*up-1);
        }
        msz/=st->GetNsites();
    }
    m_vals+=msz;
}


