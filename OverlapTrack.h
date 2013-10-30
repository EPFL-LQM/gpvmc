#ifndef _OVERLPATRACK_H
#define _OVERLAPTRACK_H
#include "Quantity.h"
#include "Stepper.h"
#include "FileManager.h"

/*!\brief class to track down the evolution of the MC weight
 *
 */
class OverlapTrack: public Quantity {
    public:
        OverlapTrack(const Stepper* stepper, FileManager* fm)
            :Quantity(stepper,fm,"OverlapTrack")
        {}
        virtual ~OverlapTrack() {}
        virtual void init()
        {
            Quantity::init();
            m_vals.clear();
        }
        virtual void save() const
        {
            double * amps=new double[m_vals.size()];
            BigDouble mean=0;
            for(size_t p=0;p<m_vals.size();++p)
                mean+=m_vals[p];
            mean/=BigDouble(double(m_vals.size()));
            for(size_t p=0;p<m_vals.size();++p)
                amps[p]=double(m_vals[p]/mean);
            m_fm->FileStream(m_basename).Write(m_vals.size(),1,amps);
            delete [] amps;
        }
        virtual void measure()
        {
            m_vals.push_back(m_stepper->weight());
        }
    protected:
        std::vector<BigDouble> m_vals;
};

#endif
