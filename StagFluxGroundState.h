#ifndef _STAGFLUXGROUNDSTATE_H
#define _STAGFLUXGROUNDSTATE_H

#include "StagFluxWaveFunction.h"
#include <vector>

class StagFluxGroundState : public StagFluxWaveFunction
{
    public:
        StagFluxGroundState(size_t Lx, size_t Ly,
                            double phi, double neel,
                            double *bc_phase)
            :StagFluxWaveFunction(Lx,Ly,
                                  Lx*Ly/2,Lx*Ly/2,
                                  phi,neel, bc_phase)
        {
            std::vector<size_t> fst(Lx*Ly,Lx*Ly/2);
            for(size_t f=0;f<Lx*Ly/2;++f)
                fst[2*f]=f;
            add_state(&fst[0],&fst[0]);
        }
};

#endif//_STAGFLUXGROUNDSTATE_H
