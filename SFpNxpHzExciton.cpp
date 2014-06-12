#include <set>
#include "SFpNxpHzExciton.h"
#include "linalg.h"

using namespace std;

SFpNxpHzExciton::SFpNxpHzExciton(size_t Lx,
                                 size_t Ly,
                                 double phi,
                                 double neel,
                                 double hz,
                                 vector<double> bc_phase,
                                 vector<size_t> q)
    :SFpNxpHzWaveFunction(Lx,Ly,
                          Lx*Ly,
                          phi,neel,hz,bc_phase),
     m_q(q)
{
    // add all two-spinons excitons
    vector<uint_vec_t> gsst(1,uint_vec_t(2*Lx*Ly,Lx*Ly));
    for(size_t f=0;f<(Lx*Ly)/2;++f){
        gsst[0][4*f]=0;
        gsst[0][4*f+1]=0;
    }
    if((m_q[0]==0 && m_q[1]==0) || (m_q[0]==Lx/2 && m_q[1]==Ly))
        AddState(gsst);
    for(size_t kx=0;kx<Lx;++kx){
        for(size_t ky=0;ky<Ly;++ky){
            if(inmbz(kx,ky)){
                size_t Qk[2]={linalg::mod(int(kx)-int(m_q[0]),m_Lx),
                              linalg::mod(int(ky)-int(m_q[1]),m_Ly)};
                mbzmod(Qk);
                for(size_t mp=0;mp<2;++mp){
                    for(size_t m=2;m<4;++m){
                        vector<uint_vec_t> st(gsst);
                        st[0][m_qn2fock[(Qk[0]*Ly+Qk[1])*4+mp]]=Lx*Ly;
                        st[0][m_qn2fock[(kx*Ly+ky)*4+m]]=0;
                        AddState(st);
                    }
                }
            }
        }
    }
}

SFpNxpHzExciton::~SFpNxpHzExciton()
{}
