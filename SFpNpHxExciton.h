#ifndef _SFPNPHXEXCITON_H
#define _SFPNPHXEXCITON_H

#include <string>
#include "SFpNpHxWaveFunction.h"

using namespace std;

class SFpNpHxExciton : public SFpNpHxWaveFunction
{
    public:
        SFpNpHxExciton(size_t Lx, size_t Ly,
                            double phi, double neel, double hx,
                            std::vector<double> bc_phase, std::vector<size_t> q);
        virtual ~SFpNpHxExciton();
    private:
        std::vector<size_t> m_q;
};

#endif//_SFPNPHXTRANSEXCITON_H

