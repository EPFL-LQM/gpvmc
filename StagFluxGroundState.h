#ifndef _STAGFLUXGROUNDSTATE_H
#define _STAGFLUXGROUNDSTATE_H

#include "StagFluxWaveFunction.h"
#include "defs.h"

class StagFluxGroundState : public StagFluxWaveFunction
{
    public:
        StagFluxGroundState(size_t Lx, size_t Ly,
                              double phi, double neel,
                              std::vector<double> bc_phase)
            :StagFluxWaveFunction(Lx,Ly,
                                    Lx*Ly/2,Lx*Ly/2,
                                    phi,neel, bc_phase)
        {
            uint_vec_t fst(Lx*Ly,Lx*Ly/2);
            for(size_t f=0;f<Lx*Ly/2;++f)
                fst[2*f]=f;
            AddState(std::vector<uint_vec_t>(2,fst));
        }
};

#endif//_STAGFLUXGROUNDSTATE_H
