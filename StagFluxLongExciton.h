#include "StagFluxWaveFunction.h"
#ifndef _STAGFLUSLONGEXCITON_H
#define _STAGFLUXLONGEXCITON_H
class StagFluxLongExciton: public StagFluxWaveFunction
{
    public:
        StagFluxLongExciton(FileManager* fm,
                            size_t Lx,size_t Ly,
                            double phi, double neel,
                            std::vector<double> bc_phase, std::vector<size_t> q);
        virtual ~StagFluxLongExciton();
    private:
        std::vector<size_t> m_q;
};
#endif//_STAGFLUXLONGEXCITON_H
