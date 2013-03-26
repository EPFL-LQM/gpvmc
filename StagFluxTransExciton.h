#ifndef _STAGFLUXTRANSEXCITON_H
#define _STAGFLUXTRANSEXCITON_H

#include <string>
#include "StagFluxWaveFunction.h"

using namespace std;

class StagFluxTransExciton : public StagFluxWaveFunction
{
    public:
        StagFluxTransExciton(size_t Lx, size_t Ly,
                             double phi, double neel,
                             double *bc_phase, size_t *q,
                             double Ecutoff=0);
        virtual ~StagFluxTransExciton();
    private:
        size_t m_q[2];
        void fock2str(const vector<size_t>& fo, string& str,size_t Nfs);
};

#endif//_STAGFLUXTRANSEXCITON_H

