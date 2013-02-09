#ifndef _MATRIXQUANTITY_H
#define _MATRIXQUANTITY_H

#include "VectorQuantity.h"

/*! \brief base class for two-dimensional (matrix) quantities.
 *
 */

class MatrixQuantity: public VectorQuantity
{
    public:
        MatrixQuantity(const Stepper* stepper,
                       FileManager* fm,
                       std::string basename,
                       int M, int N);
        std::complex<double>& Val(size_t i, size_t j)
        {
            return m_vals[i*m_N+j];
        }
        const std::complex<double>& Val(size_t i, size_t j) const
        {
            return m_vals[i*m_N+j];
        }
        virtual void save() const;
        std::string str() const;
    protected:
        size_t m_M;
        size_t m_N;
};

#endif//_MATRIXQUANTITY_H
