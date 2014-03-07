#ifndef _STAGMAGN_H
#define _STAGMAGN_H
#include "ScalarQuantity_1.h"
#include "Stepper_1.h"
#include "Amplitude_1.h"
#include "LatticeState_1.h"

class StagMagn_1: public ScalarQuantity_1 {
    public:
        StagMagn_1(const Stepper_1* stepper, FileManager* fm)
            :ScalarQuantity_1(stepper,fm,"StagMagn_1")
        {}
        virtual ~StagMagn_1()
        {}
        virtual void measure()
        {
            Quantity_1::measure();
            const LatticeState_1* st=m_stepper->GetAmp()->GetLatticeState();
            complex<double> msz(0);
            if((st->GetNfl()==1 && st->GetNifs()[0]==2) ||
                    (st->GetNfl()==2 && st->GetNifs()[0]==1 && st->GetNifs()[1]))
            {
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
                    msz+=ph*double(2*up-1);
                }
                msz/=st->GetNsites();
            }
            m_vals+=msz;
        }
};

#endif
