#ifndef _STAGFLUXGROUNDSTATE_H
#define _STAGFLUXGROUNDSTATE_H

#include "StagFluxWaveFunction.h"
#include "defs.h"

class StagFluxGroundState : public StagFluxWaveFunction
{
    public:
        StagFluxGroundState(FileManager* fm,
                            size_t Lx, size_t Ly,
                            double phi, double neel,
                            std::vector<double> bc_phase);
};

#endif//_STAGFLUXGROUNDSTATE_H
