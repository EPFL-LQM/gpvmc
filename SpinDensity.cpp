#include "SpinDensity.h"
#ifdef USEPARA
#include <omp.h>
#endif
#include "Stepper.h"
#include "Amplitude.h"
#include "SpinState.h"
#include "Jastrow.h"
#include "SpinOp.h"

using namespace std;

SpinDensity::SpinDensity(const Stepper* stepper,
                         FileManager* fm,
                         std::vector<std::vector<double> > qs,
                         size_t L)
    :MatrixQuantity(stepper,fm,"SpinDensity",2,qs.size()), m_qs(qs),
    m_ph(m_qs.size(),vector<complex<double> >(L*L,0))
    //m_xphases(m_qs.size(),std::vector<std::complex<double> >(L/2+1,0)),
    //m_yphases(m_qs.size(),std::vector<std::complex<double> >(L/2+1,0))
{
    std::complex<double> I(0,1);
    for(size_t q=0;q<m_qs.size();++q){
        for(size_t m=0;m<L;++m){
            for(size_t n=0;n<L;++n){
                m_ph[q][m*L+n]=exp(I*(m_qs[q][0]*m+m_qs[q][1]*n))/L;
                //m_xphases[q][n-m]=exp(I*m_qs[q][0]*double(int(n)-int(m)));
                //m_yphases[q][n-m]=exp(I*m_qs[q][1]*double(int(n)-int(m)));
            }
        }
    }
}

void SpinDensity::measure()
{
    const Amplitude* amp=m_stepper->GetAmp();
    const SpinState* st=amp->GetSpinState();
    Quantity::measure();
    std::vector<std::complex<double> > denum(m_qs.size(),0);
    std::vector<std::complex<double> > num(m_qs.size(),0);
    BigComplex flip(0), H_ij(0);
    size_t L=st->GetL();
    for(size_t ix=0; ix<st->GetL();++ix){
        for(size_t iy=0;iy<st->GetL();++iy){
            for(size_t lx=0;lx<st->GetL();++lx){
                for(size_t ly=0;ly<st->GetL();++ly){
                    //int tx=int(lx)-int(ix);
                    //int ty=int(ly)-int(iy);
                    //if(abs(tx)>st->GetL()/2){
                    //    if(tx>0) tx-=st->GetL();
                    //    else tx+=st->GetL();
                    //}
                    //if(abs(ty)>st->GetL()/2){
                    //    if(ty>0) ty-=st->GetL();
                    //    else ty+=st->GetL();
                    //}
                    flip=0.5*(SpinOp::SipSjm(amp,ix,iy,lx,ly)+
                              SpinOp::SimSjp(amp,ix,iy,lx,ly));
                    //bool NN=(tx==1 && ty==0) || (tx==0 && ty==1) ||
                    //        (tx==-1 && ty==0) || (tx==0 && ty==-1);
                    //if(NN){
                    //    H_ij=0.5*SpinOp::SizSjz(amp,ix,iy,lx,ly)+0.5*flip;
                    //} else {
                    //    H_ij=0.0;
                    //}
                    flip*=conj(amp->Amp());
                    //H_ij*=conj(amp->Amp());
                    if(st->GetJas()){
                        flip*=st->GetJas()->Jas();
                        //H_ij*=st->GetJas()->Jas();
                    }
                    for(size_t q=0;q<m_qs.size();++q){
                        //if(tx>0){
                        //    if(ty>0) phase=m_xphases[q][tx]*
                        //                   m_yphases[q][ty];
                        //    else phase=m_xphases[q][tx]*
                        //               conj(m_yphases[q][-ty]);
                        //} else {
                        //    if(ty>0) phase=conj(m_xphases[q][-tx])*
                        //                   m_yphases[q][ty];
                        //    else phase=conj(m_xphases[q][-tx])*
                        //               conj(m_yphases[q][-ty]);
                        //}
                        BigDouble w=m_stepper->weight();
                        denum[q]+=std::complex<double>(conj(m_ph[q][ix*L+iy])*flip*m_ph[q][lx*L+ly]/w);
                        //if(NN) num[tnum][q]+=std::complex<double>(
                        //        -4.0*(1.0-real(phase))*H_ij/w);
                    }
                }
            }
        }
    }
    for(size_t q=0;q<m_qs.size();++q){
        Val(0,q)+=num[q];
        Val(1,q)+=denum[q];
    }
}
