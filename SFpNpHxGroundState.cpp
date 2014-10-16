#include "SFpNpHxGroundState.h"

using namespace std;

SFpNpHxGroundState::SFpNpHxGroundState(FileManager* fm,
                                       size_t Lx, size_t Ly,
                                       double phi, double neel, double neel_exp, double hx,
                                       const std::vector<double>& bc_phase)
    :SFpNpHxWaveFunction(fm,Lx,Ly,
                         Lx*Ly,
                         phi,neel,neel_exp,hx,bc_phase)
{
    uint_vec_t fst(2*Lx*Ly,Lx*Ly);
    for(size_t f=0;f<Lx*Ly/2;++f){
        fst[4*f]=2*f;
        fst[4*f+1]=2*f+1;
    }
    AddState(std::vector<uint_vec_t>(1,fst));
}
