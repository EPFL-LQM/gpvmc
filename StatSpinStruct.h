#ifndef _STATSPINSTRUCT_H
#define _STATSPINSTRUCT_H
#include "MatrixQuantity.h"

class StatSpinStruct : public MatrixQuantity
{
    public:
        StatSpinStruct(const Stepper* stepper,
                       FileManager* fm);
        virtual ~StatSpinStruct() {}
        virtual void measure();
    private:
        std::vector<size_t> m_qs;
        std::vector<std::complex<double> > m_ph;
};

#endif//_STATSPINSTRUCT_H
