#ifndef _STAGMAGNZ_H
#define _STAGMAGNZ_H
#include "ScalarQuantity.h"

class StagMagnZ: public ScalarQuantity {
    public:
        StagMagnZ(const Stepper* stepper, FileManager* fm);
        virtual ~StagMagnZ();
        virtual void measure();
};

#endif
