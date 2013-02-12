#include "StagFluxWaveFunction.h"
#ifndef _STAGFLUSLONGEXCITON_H
#define _STAGFLUXLONGEXCITON_H
class StagFluxLongExciton: public StagFluxWaveFunction
{
    public:
        StagFluxLongExciton(size_t Lx, size_t Ly,
                             double phi, double neel,
                             double *bc_phase, size_t *q);
        virtual ~StagFluxLongExciton();
    private:
        size_t m_q[2];
};
#endif//_STAGFLUXLONGEXCITON_H
