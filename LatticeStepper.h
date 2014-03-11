#ifndef _SPINEXCITONSTEPPER_H
#define _SPINEXCITONSTEPPER_H

#include <algorithm>
#include "Stepper.h"
#include "BigComplex.h"
#include "BigDouble.h"
#include "defs.h"

/*! \brief Class to perform a Monte Carlo step, by swaping nearest neighbours spins.
 *
 */

class LatticeStepper: public Stepper {
    private:
        size_t m_Nfl;
        uint_vec_t m_Nifs;
        std::vector<hop_path_t> m_prev;
        int m_khop;
        BigDouble m_prev_weight;
        BigDouble m_weight;
    public:
        LatticeStepper(Amplitude* amp);
        virtual ~LatticeStepper() {}
        virtual void Reset();

        virtual BigDouble trystep();
        virtual void step();
        virtual BigDouble weight();
        virtual BigDouble weight() const {return m_weight;}
        virtual double transprob();
};

#endif//_SPINEXCITONSTEPPER_H
