#ifndef _SPINSPINCORR_H
#define _SPINSPINCORR_H
#include "MatrixQuantity.h"

/*! \brief a class to measure \f$\langle S^z(\mathbf{0})S^z(\mathbf{r})\rangle\f$
 *
 */

class SpinSpinCorr: public MatrixQuantity{
    private:
        size_t* m_R;
        size_t m_Nr;
    public:
        SpinSpinCorr(const Stepper* stepper, FileManager* fm, size_t Nr, size_t* R);
        virtual ~SpinSpinCorr();
        virtual void measure();
};


#endif//_SPINSPINCORR_H
