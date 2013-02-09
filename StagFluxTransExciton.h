#ifndef _STAGFLUXTRANSEXCITON_H
#define _STAGFLUXTRANSEXCITON_H

#include "StagFluxWaveFunction.h"

class StagFluxTransExciton : public StagFluxWaveFunction
{
    public:
        StagFluxTransExciton(size_t Lx, size_t Ly,
                             double phi, double neel,
                             double *bc_phase, size_t *q);
        virtual ~StagFluxTransExciton();
    private:
        size_t m_q[2];
};

#endif//_STAGFLUXTRANSEXCITON_H

