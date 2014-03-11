#ifndef _ScalarQuantity_H
#define _ScalarQuantity_H
#include "Quantity.h"
#include <complex>

/*! \brief Base class for scalar quantities.
 *
 */

class ScalarQuantity: public Quantity{
    public:
        ScalarQuantity(const Stepper* stepper, FileManager* fm, std::string basename)
            :Quantity(stepper,fm,basename){}
        virtual ~ScalarQuantity() {}
        virtual void init();
        virtual void save() const;
        std::string str() const;
    protected:
        std::complex<double>  m_vals;
};

#endif
