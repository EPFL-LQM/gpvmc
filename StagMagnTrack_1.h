#ifndef _STAGMAGNTRACK_1_H
#define _STAGMAGNTRACK_1_H
#include "Stepper_1.h"
#include "Amplitude_1.h"
#include "LatticeState_1.h"
#include "FileManager.h"

/*!\brief class to track down the evolution of the staggered magnetization
 *
 */
class StagMagnTrack_1: public Quantity_1 {
    public:
        StagMagnTrack_1(const Stepper_1* stepper, FileManager* fm)
            :Quantity_1(stepper,fm,"StagMagnTrack")
        {}
        virtual ~StagMagnTrack_1() {}
        virtual void init()
        {
            Quantity_1::init();
            m_vals.clear();
        }
        virtual void save() const
        {
            m_fm->FileStream(m_basename).Write(m_vals.size()*2,1,(double*)m_vals.data());
        }
        virtual void measure()
        {
            const LatticeState_1* st=m_stepper->GetAmp()->GetLatticeState();
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
    protected:
        std::vector<complex<double> > m_vals;
};

#endif
