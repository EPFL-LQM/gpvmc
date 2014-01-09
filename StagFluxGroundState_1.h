#ifndef _STAGFLUXGROUNDSTATE_1_H
#define _STAGFLUXGROUNDSTATE_1_H

#include "StagFluxWaveFunction_1.h"
#include <vector>

class StagFluxGroundState_1 : public StagFluxWaveFunction_1
{
    public:
        StagFluxGroundState_1(size_t Lx, size_t Ly,
                              double phi, double neel,
                              double *bc_phase)
            :StagFluxWaveFunction_1(Lx,Ly,
                                    Lx*Ly/2,Lx*Ly/2,
                                    phi,neel, bc_phase)
        {
            std::vector<size_t> fst(Lx*Ly,Lx*Ly/2);
            for(size_t f=0;f<Lx*Ly/2;++f)
                fst[2*f]=f;
            std::vector<const size_t*> st(2,&fst[0]);
            AddState(st);
        }
};

#endif//_STAGFLUXGROUNDSTATE_H
