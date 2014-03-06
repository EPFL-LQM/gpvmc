#ifndef _ScalarQuantity_1_H
#define _ScalarQuantity_1_H
#include "Quantity_1.h"
#include <complex>

/*! \brief Base class for scalar quantities.
 *
 */

class ScalarQuantity_1: public Quantity_1{
    public:
        ScalarQuantity_1(const Stepper_1* stepper, FileManager* fm, std::string basename)
            :Quantity_1(stepper,fm,basename){}
        virtual ~ScalarQuantity_1() {}
        virtual void init();
        virtual void save() const;
        std::string str() const;
    protected:
        std::complex<double>  m_vals;
};

#endif
