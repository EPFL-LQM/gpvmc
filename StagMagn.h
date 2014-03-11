#ifndef _STAGMAGN_H
#define _STAGMAGN_H
#include "ScalarQuantity.h"

class StagMagn: public ScalarQuantity {
    public:
        StagMagn(const Stepper* stepper, FileManager* fm);
        virtual ~StagMagn();
        virtual void measure();
};

#endif
