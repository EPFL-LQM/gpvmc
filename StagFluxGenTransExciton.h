#ifndef _STAGFLUXGENTRANSEXCITON_H
#define _STAGFLUXGENTRANSEXCITON_H

#include "StagFluxWaveFunction.h"

class StagFluxGenTransExciton : public StagFluxWaveFunction
{
    public:
        StagFluxGenTransExciton(size_t Lx, size_t Ly,
                             double phi, double neel,
                             double *bc_phase);
        virtual ~StagFluxGenTransExciton();
};

#endif//_STAGFLUXGENTRANSEXCITON_H

