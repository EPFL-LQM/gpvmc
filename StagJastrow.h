#ifndef _STAGJASTROW_H
#define _STAGJASTROW_H

#include "Jastrow.h"

class StagJastrow: public Jastrow {
    public:
        StagJastrow(SpinState* sp, double gamma)
            :Jastrow(sp),m_gamma(gamma) {init();}
        virtual ~StagJastrow() {}
    protected:
        virtual double jpot(const size_t* r);
        double m_gamma;
};

#endif//_STAGJASTROW_H
