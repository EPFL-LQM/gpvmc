#include <set>
#include "StagFluxTransExciton_1.h"
#include "linalg.h"

using namespace std;

StagFluxTransExciton_1::StagFluxTransExciton_1(size_t Lx,
                                               size_t Ly,
                                               double phi,
                                               double neel,
                                               double *bc_phase,
                                               size_t *q,
                                               double Ecutoff)
    :StagFluxWaveFunction_1(Lx,Ly,
                            Lx*Ly/2+1,Lx*Ly/2-1,
                            phi,neel,bc_phase)
{
    m_q[0]=q[0];
    m_q[1]=q[1];
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
        vector<const size_t*> st(2);
        st[0]=&stup[0];
        st[1]=&stdo[0];
        AddState(st);
    }
    // add four-spinons excitons such that the mean-field energy of
    // the excitation is less than Ecutoff
    set<string> pickedstates;
    for(size_t fe1=0;fe1<Lx*Ly/2;++fe1){
        size_t k1[2]={m_fock2qn[3*(2*fe1+1)],m_fock2qn[3*(2*fe1+1)+1]};
        for(size_t fe2=0;fe2<Lx*Ly/2;++fe2){
            size_t k2[2]={m_fock2qn[3*(2*fe2+1)],m_fock2qn[3*(2*fe2+1)+1]};
            for(size_t fe3=0;fe3<Lx*Ly/2;++fe3){
                size_t k3[2]={m_fock2qn[3*(2*fe3+1)],m_fock2qn[3*(2*fe3+1)+1]};
                size_t k4[2]={linalg::mod(int(k1[0])+int(k3[0])-int(k2[0])-int(m_q[0]),m_Lx),
                              linalg::mod(int(k1[1])+int(k3[1])-int(k2[1])-int(m_q[1]),m_Ly)};
                mbzmod(k4);
                size_t fe4=Lx*Ly/2;
                for(size_t f=0;f<Lx*Ly/2;++f){
                    if(m_fock2qn[3*(2*f+1)]==k4[0] && m_fock2qn[3*(2*f+1)+1]==k4[1]){
                        fe4=f;
                    }
                }
                if(fe4==Lx*Ly/2){
                    throw(std::logic_error("Excited state k4=k1+k3-k2-q could not be found. This is a bug."));
                }
                // calculate Delta-E. Energies associated to fe2 and fe4 are negative thus the simple sum below.
                double DE=en[fe1]+en[fe3]+en[fe2]+en[fe4];
                if(DE<Ecutoff){
                    for(size_t sig=0;sig<2;++sig){
                        // remove when either creating twice same state or destroying twice.
                        if((sig && (fe2==fe4)) || (!sig && (fe1==fe3))) continue;
                        size_t up_count=0;
                        size_t do_count=0;
                        vector<size_t> stup(Lx*Ly,Lx*Ly/2+1);
                        vector<size_t> stdo(Lx*Ly,Lx*Ly/2-1);
                        // create fock state, populating the bottom band first.
                        for(size_t f=0;f<Lx*Ly/2;++f){
                            if(sig){
                                stup[2*f]=up_count;
                                ++up_count;
                                if(f!=fe2 && f!=fe4){
                                    stdo[2*f]=do_count;
                                    ++do_count;
                                }
                            } else {
                                if(f!=fe4){
                                    stup[2*f]=up_count;
                                    ++up_count;
                                }
                                if(f!=fe2){
                                    stdo[2*f]=do_count;
                                    ++do_count;
                                }
                            }
                        }
                        // add the top band particles
                        if(sig){
                            stup[2*fe1+1]=up_count;
                            stdo[2*fe3+1]=do_count;
                        } else {
                            stup[2*fe3+1]=up_count;
                            stup[2*fe1+1]=up_count+1;
                        }
                        string strup,strdo;
                        fock2str(stup,strup,Lx*Ly/2+1);
                        fock2str(stdo,strdo,Lx*Ly/2-1);
                        if(pickedstates.find(strup+strdo)==pickedstates.end()){
                            vector<const size_t*> st(2);
                            st[0]=&stup[0];
                            st[1]=&stdo[0];
                            AddState(st);
                            pickedstates.insert(strup+strdo);
                        }
                    }
                }
            }
        }
    }
}

void StagFluxTransExciton_1::fock2str(const vector<size_t>& fo, string& str, size_t Nfs)
{
    str.resize(fo.size());
    for(size_t n=0;n<fo.size();++n){
        if(fo[n]==Nfs) str[n]='e';
        else str[n]='0';
    }
}

StagFluxTransExciton_1::~StagFluxTransExciton_1()
{}
