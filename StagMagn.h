#ifndef _STAGMAGN_H
#define _STAGMAGN_H
#include "ScalarQuantity.h"
#include "Stepper.h"
#include "Amplitude.h"
#include "SpinState.h"

class StagMagn: public ScalarQuantity {
    public:
        StagMagn(const Stepper* stepper, FileManager* fm)
            :ScalarQuantity(stepper,fm,"StagMagn")
        {}
        virtual ~StagMagn()
        {}
        virtual void measure()
        {
            Quantity::measure();
            const SpinState* st=m_stepper->GetAmp()->GetSpinState();
            double msz=0;
            for(size_t ix=0;ix<st->GetL();++ix){
                for(size_t iy=0;iy<st->GetL();++iy){
                    if(st->GetLatOc(ix,iy)==UP)
                        msz+=2*((int(ix)+int(iy))%2)-1;
                    else
                        msz-=2*((int(ix)+int(iy))%2)-1;
                }
            }
            m_vals+=msz/st->GetL()/st->GetL();
        }
};

#endif
