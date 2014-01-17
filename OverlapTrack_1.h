#ifndef _OVERLPATRACK_1_H
#define _OVERLAPTRACK_1_H
#include "Quantity_1.h"
#include "Stepper_1.h"
#include "FileManager.h"

/*!\brief class to track down the evolution of the MC weight
 *
 */
class OverlapTrack_1: public Quantity_1 {
    public:
        OverlapTrack_1(const Stepper_1* stepper, FileManager* fm)
            :Quantity_1(stepper,fm,"OverlapTrack")
        {}
        virtual ~OverlapTrack_1() {}
        virtual void init()
        {
            Quantity_1::init();
            m_vals.clear();
        }
        virtual void save() const
        {
            vector<double> amps(m_vals.size());
            BigDouble mean=0;
            for(size_t p=0;p<m_vals.size();++p)
                mean+=m_vals[p];
            mean/=BigDouble(double(m_vals.size()));
            for(size_t p=0;p<m_vals.size();++p)
                amps[p]=double(m_vals[p]/mean);
            m_fm->FileStream(m_basename).Write(m_vals.size(),1,amps.data());
        }
        virtual void measure()
        {
            m_vals.push_back(m_stepper->weight());
        }
    protected:
        std::vector<BigDouble> m_vals;
};

#endif
