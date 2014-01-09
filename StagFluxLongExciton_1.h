#include "StagFluxWaveFunction_1.h"
#ifndef _STAGFLUSLONGEXCITON_1_H
#define _STAGFLUXLONGEXCITON_1_H
class StagFluxLongExciton_1: public StagFluxWaveFunction_1
{
    public:
        StagFluxLongExciton_1(size_t L,
                              double phi, double neel,
                              double *bc_phase, size_t *q);
        virtual ~StagFluxLongExciton_1();
    private:
        size_t m_q[2];
};
#endif//_STAGFLUXLONGEXCITON_H
