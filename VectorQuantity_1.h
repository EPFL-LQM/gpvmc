#ifndef _VECTORQUANTITY_H
#define _VECTORQUANTITY_H

#include "Quantity_1.h"
#include <vector>
#include <complex>

/*! \brief Base class for vector (one-dimensional) quantities.
 *
 */

class VectorQuantity_1: public Quantity_1{
    public:
        VectorQuantity_1(const Stepper_1* stepper, FileManager* fm, std::string basename, int size);
        virtual void init();
        virtual void save() const;
        std::string str() const;
    protected:
        std::vector<std::complex<double> > m_vals;
};

#endif//_VECTORQUANTITY_H
