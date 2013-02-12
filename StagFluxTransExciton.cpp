#include "StagFluxTransExciton.h"
#include "linalg.h"

using namespace std;

StagFluxTransExciton::StagFluxTransExciton(size_t Lx,
                                           size_t Ly,
                                           double phi,
                                           double neel,
                                           double *bc_phase,
                                           size_t *q)
    :StagFluxWaveFunction(Lx,Ly,
                          Lx*Ly/2+1,Lx*Ly/2-1,
                          phi,neel,bc_phase)
{
    m_q[0]=q[0];
    m_q[1]=q[1];
    for(size_t fe=0;fe<Lx*Ly/2;++fe){
        size_t kx=m_fock2qn[3*(2*fe+1)];
        size_t ky=m_fock2qn[3*(2*fe+1)+1];
        size_t Qk[2]={linalg::mod(int(kx)-int(m_q[0]),m_Lx),
                      linalg::mod(int(ky)-int(m_q[1]),m_Ly)};
        mbzmod(Qk);
        vector<size_t> stup(Lx*Ly,Lx*Ly/2+1);
        vector<size_t> stdo(Lx*Ly,Lx*Ly/2-1);
        size_t do_count(0);
        for(size_t f=0;f<(Lx*Ly)/2;++f){
            stup[2*f]=f;
            if(m_qn2fock[(Qk[0]*Ly+Qk[1])*2]!=2*f){
                stdo[2*f]=do_count;
                ++do_count;
            }
        }
        stup[2*fe+1]=Lx*Ly/2;
        add_state(&stup[0],&stdo[0]);
    }
}

StagFluxTransExciton::~StagFluxTransExciton()
{}
