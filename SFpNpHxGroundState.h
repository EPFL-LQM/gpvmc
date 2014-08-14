#ifndef _SFPNPHXGROUNDSTATE_H
#define _SFPNPHXGROUNDSTATE_H

#include "SFpNpHxWaveFunction.h"
#include "defs.h"

class SFpNpHxGroundState : public SFpNpHxWaveFunction
{
    public:
        SFpNpHxGroundState(FileManager* fm,
                           size_t Lx, size_t Ly,
                           double phi, double neel, double hx,
                           const std::vector<double>& bc_phase);
};

#endif//_SFPNPHXGROUNDSTATE_H
