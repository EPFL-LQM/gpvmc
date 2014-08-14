#ifndef _SFPNXPHZGROUNDSTATE_H
#define _SFPNXPHZGROUNDSTATE_H

#include "SFpNxpHzWaveFunction.h"
#include "defs.h"

class SFpNxpHzGroundState : public SFpNxpHzWaveFunction
{
    public:
        SFpNxpHzGroundState(FileManager* fm,
                            size_t Lx, size_t Ly,
                            double phi, double neel, double hz,
                            const std::vector<double>& bc_phase);
};

#endif//_SFPNxPHZGROUNDSTATE_H
