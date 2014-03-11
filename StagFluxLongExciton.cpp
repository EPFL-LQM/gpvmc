#include "StagFluxLongExciton.h"
#include "linalg.h"
#include <mpi.h>

using namespace std;

StagFluxLongExciton::StagFluxLongExciton(size_t Lx,
                                             size_t Ly,
                                             double phi,
                                             double neel,
                                             vector<double> bc_phase,
                                             vector<size_t> q)
    :StagFluxWaveFunction(Lx,Ly,
                            Lx*Ly/2,Lx*Ly/2,
                            phi,neel,bc_phase),
     m_q(q)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    for(size_t fe=0;fe<Lx*Ly/2;++fe){
        size_t kx=m_fock2qn[3*(2*fe+1)];
        size_t ky=m_fock2qn[3*(2*fe+1)+1];
        size_t Qk[2]={linalg::mod(int(kx)-int(m_q[0]),m_Lx),
                      linalg::mod(int(ky)-int(m_q[1]),m_Ly)};
        mbzmod(Qk);
        vector<uint_vec_t> st(2);
        st[0]=uint_vec_t(Lx*Ly,Lx*Ly/2);
        st[1]=uint_vec_t(Lx*Ly,Lx*Ly/2);
        size_t count(0);
        for(size_t f=0;f<(Lx*Ly)/2;++f){
            st[1][2*f]=f;
            if(m_qn2fock[(Qk[0]*Ly+Qk[1])*2]!=2*f){
                st[0][2*f]=count;
                ++count;
            }
        }
        st[0][2*fe+1]=Lx*Ly/2-1;
        AddState(st);
        count=0;
        st[0]=uint_vec_t(Lx*Ly,Lx*Ly/2);
        st[1]=uint_vec_t(Lx*Ly,Lx*Ly/2);
        for(size_t f=0;f<(Lx*Ly)/2;++f){
            st[0][2*f]=f;
            if(m_qn2fock[(Qk[0]*Ly+Qk[1])*2]!=2*f){
                st[1][2*f]=count;
                ++count;
            }
        }
        st[1][2*fe+1]=Lx*Ly/2-1;
        AddState(st);
    }
    // if at gamma point or AFM Bragg peak
    // add the groundstate as it has same
    // zero momentum
    if((m_q[0]==0 && m_q[1]==0)||
       (m_q[0]==m_Lx/2 && m_q[1]==m_Ly/2)){
        vector<uint_vec_t> st(2);
        st[0]=uint_vec_t(Lx*Ly,Lx*Ly/2);
        st[1]=uint_vec_t(Lx*Ly,Lx*Ly/2);
        for(size_t f=0;f<(Lx*Ly)/2;++f){
            st[0][2*f]=f;
            st[1][2*f]=f;
        }
        AddState(st);
    }
}

StagFluxLongExciton::~StagFluxLongExciton()
{}
