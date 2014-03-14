#include <set>
#include "SFpNpHxExciton.h"
#include "linalg.h"

using namespace std;

SFpNpHxExciton::SFpNpHxExciton(size_t Lx,
                               size_t Ly,
                               double phi,
                               double neel,
                               double hx,
                               vector<double> bc_phase,
                               vector<size_t> q)
    :SFpNpHxWaveFunction(Lx,Ly,
                         Lx*Ly,
                         phi,neel,hx,bc_phase),
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
    for(size_t fe=0;fe<Lx*Ly/2;++fe){
        size_t kx=m_fock2qn[3*(2*fe+1)];
        size_t ky=m_fock2qn[3*(2*fe+1)+1];
        size_t Qk[2]={linalg::mod(int(kx)-int(m_q[0]),m_Lx),
                      linalg::mod(int(ky)-int(m_q[1]),m_Ly)};
        mbzmod(Qk);
        vector<uint_vec_t> st(gsst);
        st[0][fe*4]=Lx*Ly;
        st[0][m_qn2fock[(Qk[0]*Ly+Qk[1])*4+2]]=0;
        AddState(st);
        st=gsst;
        st[0][fe*4]=Lx*Ly;
        st[0][m_qn2fock[(Qk[0]*Ly+Qk[1])*4+3]]=0;
        AddState(st);
        st=gsst;
        st[0][fe*4+1]=Lx*Ly;
        st[0][m_qn2fock[(Qk[0]*Ly+Qk[1])*4+2]]=0;
        AddState(st);
        st=gsst;
        st[0][fe*4+1]=Lx*Ly;
        st[0][m_qn2fock[(Qk[0]*Ly+Qk[1])*4+3]]=0;
        AddState(st);
    }
}

SFpNpHxExciton::~SFpNpHxExciton()
{}
