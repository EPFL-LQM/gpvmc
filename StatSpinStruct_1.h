#ifndef _STATSPINSTRUCT_H
#define _STATSPINSTRUCT_H
#include "MatrixQuantity_1.h"
#include "defs.h"

class StatSpinStruct_1 : public MatrixQuantity_1
{
    public:
        StatSpinStruct_1(const Stepper_1* stepper,
                       FileManager* fm);
        virtual ~StatSpinStruct_1() {}
        virtual void measure();
    private:
        bool isup(const std::vector<uint_vec_t>& in);
        bool isdo(const std::vector<uint_vec_t>& in);
        std::vector<size_t> m_qs; //!< index of reciprocal space vector
        std::vector<std::complex<double> > m_ph;
};

#endif//_STATSPINSTRUCT_H
