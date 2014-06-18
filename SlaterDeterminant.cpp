#include "SlaterDeterminant.h"
#include "WaveFunction.h"
#include "LatticeState.h"
#include <limits>
#include <complex>
#include <string>
#include <sstream>
#include <iomanip>
#include "linalg.h"
#include "Timer.h"
#include "blas_lapack.h"
#ifdef USEMPI
#include <mpi.h>
#endif

using namespace std;

#define EPS -8
#define DEPS 1e-8

SlaterDeterminant::SlaterDeterminant(LatticeState* sp, WaveFunction* wav)
    :m_lst(sp), m_wav(wav),
    m_amp(0), m_amp_ok(false)
{
    bool compat=true;
    if(m_lst->GetNfl()!=m_wav->GetNfl()){
        compat=false;
    } else {
        if(!equal(m_lst->GetNpt().begin(),m_lst->GetNpt().end(),m_wav->GetNpt().begin())){
            compat=false;
        }
    }
    if(!compat){
#ifdef EXCEPT
        throw(std::logic_error("SlaterDeterminant::"
                    "SlaterDeterminant(const LatticeState*, const WaveFunction*):"
                    "cannot define amplitude with different number of "
                    "flavours or number of particles."));
#else
        cerr<<"SlaterDeterminant::"
              "SlaterDeterminant(const LatticeState*, const WaveFunction*):"
              "cannot define amplitude with different number of "
              "flavours or number of particles."<<endl;
        abort();
#endif
    }
    m_Nfl=m_lst->GetNfl();
    m_N=m_lst->GetNpt();
    m_mat.resize(m_Nfl);
    m_mati.resize(m_Nfl);
    for(size_t fl=0;fl<m_Nfl;++fl){
        m_mat[fl].resize(m_N[fl]*m_N[fl]);
        m_mati[fl].resize(m_N[fl]*m_N[fl]);
    }
    Init();
}

SlaterDeterminant::~SlaterDeterminant()
{}

void SlaterDeterminant::Init()
{
#ifdef PROFILE
    Timer::tic("SlaterDeterminant::Init");
#endif
#ifdef DEBUG
    std::cout<<"SlaterDeterminant::Init: has been called"<<std::endl;
#endif
    for(size_t fl=0;fl<m_Nfl;++fl){
        for(size_t r=0;r<m_N[fl];++r){
            for(size_t f=0;f<m_N[fl];++f){
                m_mat[fl][r*m_N[fl]+f]=m_wav->MatEl(m_wav->Getpt()[fl][f],m_lst->Getpt()[fl][r],fl);
            }
        }
    }
    vector<BigComplex> d(m_Nfl,0.0);
    for(size_t fl=0;fl<m_Nfl;++fl)
        linalg::DetInv(m_mat[fl].data(),m_mati[fl].data(),m_N[fl],d[fl]);
    if(any_of(d.begin(),d.end(),[](const BigComplex& n){return (n.exp()<EPS) || (n==0.0);})){
        m_amp_ok=false;
        m_amp=0.0;
    } else {
        // real space state sign omitted since it is canceled
        // in this calculation with the sign of the measured
        // quantity matrix elements
        m_amp=d[0];
        for(size_t fl=1;fl<m_Nfl;++fl) m_amp*=d[fl];
        m_amp_ok=true;
    }
#ifdef PROFILE
    Timer::toc("SlaterDeterminant::Init");
#endif
}

BigComplex SlaterDeterminant::Amp() const {return m_amp*m_wav->GetSign();}

void SlaterDeterminant::VirtUpdate(const vector<vector<hop_path_t> >& rhop,
                             const vector<vector<hop_path_t> >& khop,
                             vector<BigComplex>& qs) const
{
    if(!m_amp_ok){
#ifdef EXCEPT
        throw(std::logic_error("SlaterDeterminant::VirtUpdate:"
                               " Should not be called from"
                               " a state without overlap."));
#else
        cerr<<"SlaterDeterminant::VirtUpdate:"
              " Should not be called from"
              " a state without overlap."<<endl;
        abort();
#endif
    }
    size_t Nr=rhop.size();
    size_t Nk=khop.size();
    if(Nr*Nk==0){
#ifdef EXCEPT
        throw(std::logic_error("SlaterDeterminant::VirtUpdate:"
                               " the condition min(Nr)=1 "
                               " and min(Nk)=1"
                               " must be fullfilled."));
#else
        cerr<<"SlaterDeterminant::VirtUpdate:"
              " the condition min(Nr)=1 "
              " and min(Nk)=1"
              " must be fullfilled."<<endl;
        abort();
    }
#endif
#ifdef PROFILE
    Timer::tic("SlaterDeterminant::VirtUpdate");
#endif
    size_t NN=max(Nr,Nk);
    qs.resize(Nr*Nk);
    for(size_t nk=0;nk<Nk;++nk){
        int sign=m_wav->HopSign(khop[nk]);
        for(size_t nr=0;nr<Nr;++nr){
            qs[nk*Nr+nr]=m_amp*m_wav->GetSign()*sign;
        }
    }
    vector<BigComplex> qq(NN);
    for(size_t fl=0;fl<m_Nfl;++fl){
        vector<size_t> ridx(Nr,0);
        vector<size_t> kidx(Nk,0);
        for(size_t n=1;n<Nr;++n)
            ridx[n]=ridx[n-1]+rhop[n-1][fl].size();
        for(size_t n=1;n<Nk;++n)
            kidx[n]=kidx[n-1]+khop[n-1][fl].size();
        size_t Nrf=ridx.back()+rhop.back()[fl].size();
        size_t Nkf=kidx.back()+khop.back()[fl].size();
        // allocate block matrices of change.
        vector<complex<double> > V(m_N[fl]*Nrf);
        vector<complex<double> > U(m_N[fl]*Nkf);
        /* copy column and row changes as they would
         * happen separatly (crossings between row
         * and column changes are treated later
         * below).
         */
        for(size_t n=0;n<Nr;++n){
            for(size_t r=0;r<rhop[n][fl].size();++r){
                for(size_t f=0;f<m_N[fl];++f)
                    V[(ridx[n]+r)*m_N[fl]+f]=
                        m_wav->MatEl(m_wav->Getpt()[fl][f],
                                     rhop[n][fl][r].second,fl);
            }
        }
        for(size_t n=0;n<Nk;++n){
            for(size_t r=0;r<khop[n][fl].size();++r){
                for(size_t f=0;f<m_N[fl];++f){
                    U[f*Nkf+kidx[n]+r]=
                        m_wav->MatEl(khop[n][fl][r].second,
                                     m_lst->Getpt()[fl][f],fl);
                }
            }
        }
        vector<size_t> rfid(Nrf);
        vector<size_t> kfid(Nkf);
        vector<size_t> rfr(Nr);
        vector<size_t> kfr(Nk);
        for(size_t nr=0;nr<Nr;++nr){
            rfr[nr]=rhop[nr][fl].size();
            for(size_t rr=0;rr<rhop[nr][fl].size();++rr)
                rfid[ridx[nr]+rr]=m_lst->Getfs()[fl][rhop[nr][fl][rr].first];
        }
        for(size_t nk=0;nk<Nk;++nk){
            kfr[nk]=khop[nk][fl].size();
            for(size_t rk=0;rk<khop[nk][fl].size();++rk){
                kfid[kidx[nk]+rk]=m_wav->Getfs()[fl][khop[nk][fl][rk].first];
            }
        }
        /* treat row and columns crossings and then update
         * the determinants. As linalg::DetUpdate is efficient
         * for large multiple updates, I use NN=max(Nk,Nr) for
         * simultaneous update.
         */
        if(Nk==NN){
            for(size_t nr=0;nr<Nr;++nr){
                for(size_t rr=0;rr<rhop[nr][fl].size();++rr){
                    for(size_t nk=0;nk<Nk;++nk){
                        for(size_t rk=0;rk<khop[nk][fl].size();++rk){
                            U[rfid[ridx[nr]+rr]*Nkf+kidx[nk]+rk]=
                                m_wav->MatEl(khop[nk][fl][rk].second,
                                             rhop[nr][fl][rr].second,fl);
                        }
                    }
                }
#ifdef PROFILE
                Timer::tic("SlaterDeterminant::VirtUpdate/DetUpdate");
#endif
                linalg::DetUpdate(m_mat[fl].data(),m_mati[fl].data(),m_N[fl],
                                  &V[ridx[nr]*m_N[fl]],m_N[fl],
                                  &rfid[ridx[nr]],&rfr[nr],1,
                                  U.data(),Nkf,&kfid[0],&kfr[0],Nk,
                                  qq.data());
#ifdef PROFILE
                Timer::toc("SlaterDeterminant::VirtUpdate/DetUpdate");
#endif
                for(size_t nk=0;nk<Nk;++nk){
                    qs[nk*Nr+nr]*=qq[nk];
                }
                for(size_t rr=0;rr<rhop[nr][fl].size();++rr){
                    for(size_t nk=0;nk<Nk;++nk){
                        for(size_t rk=0;rk<khop[nk][fl].size();++rk){
                            U[rfid[ridx[nr]+rr]*Nkf+kidx[nk]+rk]=
                                m_wav->MatEl(khop[nk][fl][rk].second,
                                         m_lst->Getpt()[fl][rfid[ridx[nr]+rr]],fl);
                        }
                    }
                }
            }
        } else {
            for(size_t nk=0;nk<Nk;++nk){
                for(size_t rk=0;rk<khop[nk][fl].size();++rk){
                    for(size_t nr=0;nr<Nr;++nr){
                        for(size_t rr=0;rr<rhop[nr][fl].size();++rr){
                            V[(ridx[nr]+rr)*m_N[fl]+kfid[kidx[nk]+rk]]=
                                m_wav->MatEl(khop[nk][fl][rk].second,
                                             rhop[nr][fl][rr].second,fl);
                        }
                    }
                }
#ifdef PROFILE
                Timer::tic("SlaterDeterminant::VirtUpdate/DetUpdate");
#endif
                linalg::DetUpdate(m_mat[fl].data(),m_mati[fl].data(),m_N[fl],
                                  V.data(),m_N[fl],&rfid[0],&rfr[0],Nr,
                                  &U[kidx[nk]],Nkf,&kfid[kidx[nk]],
                                  &kfr[nk],1,qq.data());
#ifdef PROFILE
                Timer::toc("SlaterDeterminant::VirtUpdate/DetUpdate");
#endif
                for(size_t nr=0;nr<Nr;++nr){
                    qs[nk*Nr+nr]*=qq[nr];
                }
                for(size_t rk=0;rk<khop[nk][fl].size();++rk){
                    for(size_t nr=0;nr<Nr;++nr){
                        for(size_t rr=0;rr<rhop[nr][fl].size();++rr){
                            V[(ridx[nr]+rr)*m_N[fl]+kfid[kidx[nk]+rk]]=
                                m_wav->MatEl(m_wav->Getpt()[fl][kfid[kidx[nk]+rk]],
                                             rhop[nr][fl][rr].second,fl);
                        }
                    }
                }
            }
        }
    }
#ifdef PROFILE
    Timer::toc("SlaterDeterminant::VirtUpdate");
#endif
}

void SlaterDeterminant::Update(const vector<hop_path_t>& rhop,
                         const vector<hop_path_t>& khop)
{
    if(!m_amp_ok){
#ifdef EXCEPT
        throw(std::logic_error("SlaterDeterminant::Update:"
                               " Should not be called from"
                               " a state without overlap."));
#else
        cerr<<"SlaterDeterminant::Update:"
              " Should not be called from"
              " a state without overlap."<<endl;
        abort();
#endif
    }
#ifdef PROFILE
    Timer::tic("SlaterDeterminant::Update");
#endif
    for(size_t fl=0;fl<m_Nfl;++fl){
        vector<complex<double> > V(m_N[fl]*rhop[fl].size());
        vector<complex<double> > U(m_N[fl]*khop[fl].size());
        for(size_t r=0;r<rhop[fl].size();++r){
            for(size_t f=0;f<m_N[fl];++f){
                V[r*m_N[fl]+f]=m_wav->MatEl(m_wav->Getpt()[fl][f],
                                            rhop[fl][r].second,fl);
            }
        }
        for(size_t r=0;r<khop[fl].size();++r){
            for(size_t f=0;f<m_N[fl];++f){
                U[f*khop[fl].size()+r]=m_wav->MatEl(khop[fl][r].second,
                                                    m_lst->Getpt()[fl][f],fl);
            }
        }
        BigComplex qq;
        vector<size_t> kid(khop[fl].size());
        vector<size_t> rid(rhop[fl].size());
        for(size_t r=0;r<khop[fl].size();++r){
            kid[r]=m_wav->Getfs()[fl][khop[fl][r].second];
        }
        for(size_t r=0;r<rhop[fl].size();++r){
            rid[r]=m_lst->Getfs()[fl][rhop[fl][r].second];
        }
#ifdef PROFILE
        Timer::tic("SlaterDeterminant::Update/InvUpdate");
#endif
        linalg::InvUpdate(m_mat[fl].data(),m_mati[fl].data(),m_N[fl],
                          V.data(),rid.data(),rhop[fl].size(),
                          U.data(),kid.data(),khop[fl].size(),qq);
#ifdef PROFILE
        Timer::toc("SlaterDeterminant::Update/InvUpdate");
#endif
        m_amp*=qq;
    }
#ifdef PROFILE
    Timer::toc("SlaterDeterminant::Update");
#endif
}

