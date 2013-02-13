#include "StagFluxGenTransExciton.h"
#include "linalg.h"

using namespace std;

StagFluxGenTransExciton::StagFluxGenTransExciton(size_t Lx,
                                                 size_t Ly,
                                                 double phi,
                                                 double neel,
                                                 double *bc_phase)
    :StagFluxWaveFunction(Lx,Ly,
                          Lx*Ly/2+1,Lx*Ly/2-1,
                          phi,neel,bc_phase)
{
    for(size_t fe1=0;fe1<Lx*Ly/2;++fe1){
        for(size_t fe2=0;fe2<Lx*Ly/2;++fe2){
            vector<size_t> stup(Lx*Ly,Lx*Ly/2+1);
            vector<size_t> stdo(Lx*Ly,Lx*Ly/2-1);
            size_t do_count(0);
            for(size_t f=0;f<(Lx*Ly)/2;++f){
                stup[2*f]=f;
                if(f!=fe2){
                    stdo[2*f]=do_count;
                    ++do_count;
                }
            }
            stup[2*fe1+1]=Lx*Ly/2;
            add_state(&stup[0],&stdo[0]);
        }
    }
}

StagFluxGenTransExciton::~StagFluxGenTransExciton()
{}
