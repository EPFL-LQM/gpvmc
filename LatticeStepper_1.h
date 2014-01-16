#ifndef _SPINEXCITONSTEPPER_1_H
#define _SPINEXCITONSTEPPER_1_H

#include <algorithm>
#include "linalg.h"
#include "RanGen.h"
#include "Timer.h"
#include "Stepper_1.h"
#include "BigComplex.h"
#include "BigDouble.h"
#include "defs.h"

/*! \brief Class to perform a Monte Carlo step, by swaping nearest neighbours spins.
 *
 */

class LatticeStepper_1: public Stepper_1 {
    private:
        size_t m_Nfl;
        uint_vec_t m_Nifs;
        std::vector<hop_path_t> m_prev;
        int m_khop;
        BigDouble m_prev_weight;
        BigDouble m_weight;
    public:
        LatticeStepper_1(Amplitude_1* amp);
        virtual ~LatticeStepper_1() {}
        virtual void Reset();

        virtual BigDouble trystep();
        virtual void step();
        virtual BigDouble weight();
        virtual BigDouble weight() const {return m_weight;}
        virtual double transprob();
};

#endif//_SPINEXCITONSTEPPER_H
