#ifndef _STAGFLUXTRANSEXCITON_H
#define _STAGFLUXTRANSEXCITON_H

#include <string>
#include "StagFluxWaveFunction.h"

using namespace std;

class StagFluxTransExciton : public StagFluxWaveFunction
{
    public:
        StagFluxTransExciton(FileManager* fm,
                             size_t Lx, size_t Ly,
                             double phi, double neel,
                             std::vector<double> bc_phase, std::vector<size_t> q);
        virtual ~StagFluxTransExciton();
    private:
        std::vector<size_t> m_q;
        void fock2str(const vector<size_t>& fo, string& str,size_t Nfs);
};

#endif//_STAGFLUXTRANSEXCITON_H

