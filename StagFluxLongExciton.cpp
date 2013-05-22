#include "StagFluxLongExciton.h"
#include "linalg.h"
#include <mpi.h>

using namespace std;

StagFluxLongExciton::StagFluxLongExciton(size_t L,
                                         double phi,
                                         double neel,
                                         double *bc_phase,
                                         size_t *q)
    :StagFluxWaveFunction(L,L,
                          L*L/2,L*L/2,
                          phi,neel,bc_phase)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    m_q[0]=q[0];
    m_q[1]=q[1];
    for(size_t fe=0;fe<L*L/2;++fe){
        size_t kx=m_fock2qn[3*(2*fe+1)];
        size_t ky=m_fock2qn[3*(2*fe+1)+1];
        size_t Qk[2]={linalg::mod(int(kx)-int(m_q[0]),m_Lx),
                      linalg::mod(int(ky)-int(m_q[1]),m_Ly)};
        mbzmod(Qk);
        vector<size_t> stup(L*L,L*L/2);
        vector<size_t> stdo(L*L,L*L/2);
        size_t count(0);
        for(size_t f=0;f<(L*L)/2;++f){
            stdo[2*f]=f;
            if(m_qn2fock[(Qk[0]*L+Qk[1])*2]!=2*f){
                stup[2*f]=count;
                ++count;
            }
        }
        stup[2*fe+1]=L*L/2-1;
        add_state(&stup[0],&stdo[0]);
        count=0;
        stup=vector<size_t>(L*L,L*L/2);
        stdo=vector<size_t>(L*L,L*L/2);
        for(size_t f=0;f<(L*L)/2;++f){
            stup[2*f]=f;
            if(m_qn2fock[(Qk[0]*L+Qk[1])*2]!=2*f){
                stdo[2*f]=count;
                ++count;
            }
        }
        stdo[2*fe+1]=L*L/2-1;
        add_state(&stup[0],&stdo[0]);
    }
    // if at gamma point or AFM Bragg peak
    // add the groundstate as it has same
    // zero momentum
    if((m_q[0]==0 && m_q[1]==0)||
       (m_q[0]==m_Lx/2 && m_q[1]==m_Ly/2)){
        vector<size_t> stup(L*L,L*L/2);
        vector<size_t> stdo(L*L,L*L/2);
        for(size_t f=0;f<(L*L)/2;++f){
            stup[2*f]=f;
            stdo[2*f]=f;
        }
        add_state(&stup[0],&stdo[0]);
    }
}

StagFluxLongExciton::~StagFluxLongExciton()
{}
