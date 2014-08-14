#include <set>
#include "StagFluxTransExciton.h"
#include "linalg.h"

using namespace std;

StagFluxTransExciton::StagFluxTransExciton(FileManager* fm,
                                           size_t Lx,
                                           size_t Ly,
                                           double phi,
                                           double neel,
                                           vector<double> bc_phase,
                                           vector<size_t> q)
    :StagFluxWaveFunction(fm,Lx,Ly,
                          Lx*Ly/2+1,Lx*Ly/2-1,
                          phi,neel,bc_phase),
     m_q(q)
{
    vector<double> en(Lx*Ly/2,0);
    // add all two-spinons excitons
    for(size_t fe=0;fe<Lx*Ly/2;++fe){
        size_t kx=m_fock2qn[3*(2*fe+1)];
        size_t ky=m_fock2qn[3*(2*fe+1)+1];
        double k[2]={double(kx)*2*M_PI/Lx,double(ky)*2*M_PI/Ly};
        en[fe]=omega(k);
        size_t Qk[2]={linalg::mod(int(kx)-int(m_q[0]),m_Lx),
                      linalg::mod(int(ky)-int(m_q[1]),m_Ly)};
        mbzmod(Qk);
        vector<uint_vec_t> st(2);
        st[0]=uint_vec_t(Lx*Ly,Lx*Ly/2+1);
        st[1]=uint_vec_t(Lx*Ly,Lx*Ly/2-1);
        size_t do_count(0);
        for(size_t f=0;f<(Lx*Ly)/2;++f){
            st[0][2*f]=f;
            if(m_qn2fock[(Qk[0]*Ly+Qk[1])*2]!=2*f){
                st[1][2*f]=do_count;
                ++do_count;
            }
        }
        st[0][2*fe+1]=Lx*Ly/2;
        AddState(st);
    }
}

void StagFluxTransExciton::fock2str(const vector<size_t>& fo, string& str, size_t Nfs)
{
    str.resize(fo.size());
    for(size_t n=0;n<fo.size();++n){
        if(fo[n]==Nfs) str[n]='e';
        else str[n]='0';
    }
}

StagFluxTransExciton::~StagFluxTransExciton()
{}
