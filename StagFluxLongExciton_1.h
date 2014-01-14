#include "StagFluxWaveFunction_1.h"
#ifndef _STAGFLUSLONGEXCITON_1_H
#define _STAGFLUXLONGEXCITON_1_H
class StagFluxLongExciton_1: public StagFluxWaveFunction_1
{
    public:
        StagFluxLongExciton_1(size_t Lx,size_t Ly,
                              double phi, double neel,
                              std::vector<double> bc_phase, std::vector<size_t> q);
        virtual ~StagFluxLongExciton_1();
    private:
        std::vector<size_t> m_q;
};
#endif//_STAGFLUXLONGEXCITON_H
