#ifndef _STAGMAGN_H
#define _STAGMAGN_H
#include "VectorQuantity.h"

class StagMagn: public VectorQuantity {
    public:
        StagMagn(const Stepper* stepper, FileManager* fm, bool meas_trans=false);
        virtual ~StagMagn();
        virtual void measure();
    private:
        bool m_meas_trans;
};

#endif
