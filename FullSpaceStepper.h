#ifndef _SPINEXCITONSTEPPER_H
#define _SPINEXCITONSTEPPER_H

#include <algorithm>
#include "linalg.h"
#include "RanGen.h"
#include "Timer.h"
#include "Stepper.h"
#include "BigComplex.h"
#include "BigDouble.h"

/*! \brief Class to perform a Monte Carlo step, by swaping nearest neighbours spins.
 *
 */

class FullSpaceStepper: public Stepper {
    private:
        int m_prev_up;
        int m_prev_down;
        int m_khop;
        BigDouble m_prev;
        BigDouble m_weight;
    public:
        FullSpaceStepper(Amplitude* amp)
            :Stepper(amp), m_prev_up(-1), m_prev_down(-1), m_khop(-1), m_weight(-1){}
        virtual ~FullSpaceStepper() {}
        virtual void Reset();

        virtual BigDouble trystep();
        virtual void step();
        virtual BigDouble weight();
        virtual BigDouble weight() const {return m_weight;}
        virtual double transprob();
};

#endif//_SPINEXCITONSTEPPER_H
