#ifndef _STEPPER_1_H
#define _STEPPER_1_H
#include "BigDouble.h"

class Amplitude_1;

/*! \brief Class to perform a Monte Carlo step, by swaping nearest neighbours spins.
 *
 */

class Stepper_1 {
    public:
        Stepper_1(Amplitude_1* amp)
            : m_amp(amp){}
        virtual ~Stepper_1() {}
        virtual void Reset()=0;

        virtual BigDouble trystep()=0;
        virtual void step()=0;
        virtual BigDouble weight()=0;
        virtual BigDouble weight() const=0;
        virtual double transprob()=0;

        const Amplitude_1* GetAmp() const {return m_amp;}
        Amplitude_1* GetAmp() {return m_amp;}

    protected:
        Amplitude_1* m_amp;
};

#endif//_STEPPER_H
