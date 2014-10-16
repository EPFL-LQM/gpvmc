#ifndef _SFPNXPHZEXCITON_H
#define _SFPNXPHZEXCITON_H

#include <string>
#include "SFpNxpHzWaveFunction.h"

using namespace std;

class SFpNxpHzExciton : public SFpNxpHzWaveFunction
{
    public:
        SFpNxpHzExciton(FileManager* fm,
                        size_t Lx, size_t Ly,
                        double phi, double neel,
                        double neel_exp, double hz,
                        std::vector<double> bc_phase, std::vector<size_t> q);
        virtual ~SFpNxpHzExciton();
    private:
        std::vector<size_t> m_q;
};

#endif//_SFPNXPHZEXCITON_H

