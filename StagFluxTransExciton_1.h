#ifndef _STAGFLUXTRANSEXCITON_1_H
#define _STAGFLUXTRANSEXCITON_1_H

#include <string>
#include "StagFluxWaveFunction_1.h"

using namespace std;

class StagFluxTransExciton_1 : public StagFluxWaveFunction_1
{
    public:
        StagFluxTransExciton_1(size_t Lx, size_t Ly,
                               double phi, double neel,
                               double *bc_phase, size_t *q,
                               double Ecutoff=0);
        virtual ~StagFluxTransExciton_1();
    private:
        size_t m_q[2];
        void fock2str(const vector<size_t>& fo, string& str,size_t Nfs);
};

#endif//_STAGFLUXTRANSEXCITON_H

