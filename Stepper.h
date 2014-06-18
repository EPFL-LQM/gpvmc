#ifndef _STEPPER_H
#define _STEPPER_H
#include "BigDouble.h"

class SlaterDeterminant;
class Jastrow;

/*! \brief Class to perform a Monte Carlo step, by swaping nearest neighbours spins.
 *
 */

class Stepper {
    public:
        Stepper(SlaterDeterminant* amp, Jastrow* jas)
            : m_amp(amp),m_jas(jas){}
        virtual ~Stepper() {}
        virtual void Reset()=0;

        virtual BigDouble trystep()=0;
        virtual void step()=0;
        virtual BigDouble weight()=0;
        virtual BigDouble weight() const=0;
        virtual double transprob()=0;

        const SlaterDeterminant* GetAmp() const {return m_amp;}
        SlaterDeterminant* GetAmp() {return m_amp;}
        const Jastrow* GetJas() const {return m_jas;}
        Jastrow* GetJas() {return m_jas;}

    protected:
        SlaterDeterminant* m_amp;
        Jastrow* m_jas;
};

#endif//_STEPPER_H
