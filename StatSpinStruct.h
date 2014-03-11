#ifndef _STATSPINSTRUCT_H
#define _STATSPINSTRUCT_H
#include "MatrixQuantity.h"
#include "defs.h"

class StatSpinStruct : public MatrixQuantity
{
    public:
        StatSpinStruct(const Stepper* stepper,
                       FileManager* fm);
        virtual ~StatSpinStruct() {}
        virtual void measure();
    private:
        bool isup(const std::vector<uint_vec_t>& in);
        bool isdo(const std::vector<uint_vec_t>& in);
        std::vector<size_t> m_qs; //!< index of reciprocal space vector
        std::vector<std::complex<double> > m_ph;
};

#endif//_STATSPINSTRUCT_H
