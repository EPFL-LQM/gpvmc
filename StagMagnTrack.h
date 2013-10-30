#ifndef _STAGMAGNTRACK_H
#define _STAGMAGNTRACK_H
#include "ScalarQuantity.h"
#include "Stepper.h"
#include "Amplitude.h"
#include "SpinState.h"
#include "FileManager.h"

/*!\brief class to track down the evolution of the staggered magnetization
 *
 */
class StagMagnTrack: public Quantity {
    public:
        StagMagnTrack(const Stepper* stepper, FileManager* fm)
            :Quantity(stepper,fm,"StagMagnTrack")
        {}
        virtual ~StagMagnTrack() {}
        virtual void init()
        {
            Quantity::init();
            m_vals.clear();
        }
        virtual void save() const
        {
            m_fm->FileStream(m_basename).Write(m_vals.size(),1,&m_vals[0]);
        }
        virtual void measure()
        {
            m_vals.push_back(0);
            const SpinState* st=m_stepper->GetAmp()->GetSpinState();
            for(size_t x=0;x<st->GetL();++x){
                for(size_t y=0;y<st->GetL();++y){
                    int phase=1-2*(int(x+y)%2);
                    if(st->GetLatOc(x,y)==UP)
                        m_vals.back()+=phase;
                    else if(st->GetLatOc(x,y)==DOWN)
                        m_vals.back()-=phase;
                }
            }
            m_vals.back()/=pow(st->GetL(),2);
        }
    protected:
        std::vector<double> m_vals;
};

#endif
