#ifndef _STEPPER_H
#define _STEPPER_H
#include "BigDouble.h"

class Amplitude;

/*! \brief Class to perform a Monte Carlo step, by swaping nearest neighbours spins.
 *
 */

class Stepper {
    public:
        Stepper(Amplitude* amp)
            : m_amp(amp){}
        virtual ~Stepper() {}
        virtual void Reset()=0;

        virtual BigDouble trystep()=0;
        virtual void step()=0;
        virtual BigDouble weight()=0;
        virtual BigDouble weight() const=0;
        virtual double transprob()=0;

        const Amplitude* GetAmp() const {return m_amp;}
        Amplitude* GetAmp() {return m_amp;}

    protected:
        Amplitude* m_amp;
};

#endif//_STEPPER_H
