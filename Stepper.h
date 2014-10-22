#ifndef _STEPPER_H
#define _STEPPER_H
#include "BigDouble.h"

class SlaterDeterminant;
class LatticeState;
class WaveFunction;

/*! \brief Class to perform a Monte Carlo step, by swaping nearest neighbours spins.
 *
 */

class Stepper {
    public:
        Stepper(LatticeState* latstate, WaveFunction* wav, SlaterDeterminant* amp)
            : m_amp(amp),m_latstate(latstate),m_wav(wav) {}
        virtual ~Stepper() {}
        virtual void Reset()=0;

        virtual BigDouble trystep()=0;
        virtual void step()=0;
        virtual BigDouble weight()=0;
        virtual BigDouble weight() const=0;
        virtual double transprob()=0;

        const SlaterDeterminant* GetAmp() const {return m_amp;}
        const LatticeState* GetLatticeState() const {return m_latstate;}
        const WaveFunction* GetWaveFunction() const {return m_wav;}

    protected:
        SlaterDeterminant* m_amp;
        LatticeState* m_latstate;
        WaveFunction* m_wav;
};

#endif//_STEPPER_H
