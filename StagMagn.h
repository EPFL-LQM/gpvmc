#ifndef _STAGMAGN_H
#define _STAGMAGN_H
#include "VectorQuantity.h"

class StagMagn: public VectorQuantity {
    public:
        StagMagn(const Stepper* stepper, FileManager* fm);
        virtual ~StagMagn();
        virtual void measure();
};

#endif
