#include "SpinDensity.h"
#ifdef USEPARA
#include <omp.h>
#endif
#include "Stepper.h"
#include "Amplitude.h"
#include "SpinState.h"
#include "Jastrow.h"
#include "SpinOp.h"

SpinDensity::SpinDensity(const Stepper* stepper,
                         FileManager* fm,
                         std::vector<std::vector<double> > qs,
                         size_t L)
    :MatrixQuantity(stepper,fm,"SpinDensity",2,qs.size()), m_qs(qs),
    m_xphases(m_qs.size(),std::vector<std::complex<double> >(L/2+1,0)),
    m_yphases(m_qs.size(),std::vector<std::complex<double> >(L/2+1,0))
{
    std::complex<double> I(0,1);
    for(size_t q=0;q<m_qs.size();++q){
        for(size_t m=0;m<=L/2;++m){
            for(size_t n=m;n<=L/2;++n){
                m_xphases[q][n-m]=exp(I*m_qs[q][0]*double(int(n)-int(m)));
                m_yphases[q][n-m]=exp(I*m_qs[q][1]*double(int(n)-int(m)));
            }
        }
    }
}

void SpinDensity::measure()
{
    const Amplitude* amp=m_stepper->GetAmp();
    const SpinState* st=amp->GetSpinState();
    Quantity::measure();
    std::vector<std::vector<std::complex<double> > > denum;
    std::vector<std::vector<std::complex<double> > > num;
    int num_threads=1;
    size_t tnum=0;
#ifdef USEPARA
#pragma omp parallel
    {
        num_threads=omp_get_num_threads();
        tnum=omp_get_thread_num();
#pragma omp single
        {
#endif
            denum=std::vector<std::vector<std::complex<double> > >(num_threads,
                    std::vector<std::complex<double> >(m_qs.size(),0));
            num=std::vector<std::vector<std::complex<double> > >(num_threads,
                    std::vector<std::complex<double> >(m_qs.size(),0));
#ifdef USEPARA
        }
#pragma omp barrier
#endif
        std::complex<double> phase;
        BigComplex flip(0), H_ij(0);
#ifdef USEPARA
        size_t chunk=st->GetL()/omp_get_num_threads();
        if(omp_get_num_threads()*chunk<st->GetL())++chunk;
#pragma omp for schedule(static,chunk)
#endif
        for(size_t ix=0; ix<st->GetL();++ix){
            for(size_t iy=0;iy<st->GetL();++iy){
                for(size_t lx=0;lx<st->GetL();++lx){
                    for(size_t ly=0;ly<st->GetL();++ly){
                        int tx=int(lx)-int(ix);
                        int ty=int(ly)-int(iy);
                        if(abs(tx)>st->GetL()/2){
                            if(tx>0) tx-=st->GetL();
                            else tx+=st->GetL();
                        }
                        if(abs(ty)>st->GetL()/2){
                            if(ty>0) ty-=st->GetL();
                            else ty+=st->GetL();
                        }
                        flip=0.5*(SpinOp::SipSjm(amp,ix,iy,lx,ly)+
                                  SpinOp::SimSjp(amp,ix,iy,lx,ly));
                        bool NN=(tx==1 && ty==0) || (tx==0 && ty==1) ||
                                (tx==-1 && ty==0) || (tx==0 && ty==-1);
                        if(NN){
                            H_ij=0.5*SpinOp::SizSjz(amp,ix,iy,lx,ly)+0.5*flip;
                        } else {
                            H_ij=0.0;
                        }
                        flip*=conj(amp->Amp());
                        H_ij*=conj(amp->Amp());
                        if(st->GetJas()){
                            flip*=st->GetJas()->Jas();
                            H_ij*=st->GetJas()->Jas();
                        }
                        for(size_t q=0;q<m_qs.size();++q){
                            if(tx>0){
                                if(ty>0) phase=m_xphases[q][tx]*
                                               m_yphases[q][ty];
                                else phase=m_xphases[q][tx]*
                                           conj(m_yphases[q][-ty]);
                            } else {
                                if(ty>0) phase=conj(m_xphases[q][-tx])*
                                               m_yphases[q][ty];
                                else phase=conj(m_xphases[q][-tx])*
                                           conj(m_yphases[q][-ty]);
                            }
                            BigDouble w=m_stepper->weight();
                            denum[tnum][q]+=std::complex<double>(phase*flip/w);
                            if(NN) num[tnum][q]+=std::complex<double>(
                                    -4.0*(1.0-real(phase))*H_ij/w);
                        }
                    }
                }
            }
        }
#ifdef USEPARA
    }
#endif
    for(size_t t=0;t<num.size();++t){
        for(size_t q=0;q<m_qs.size();++q){
            Val(0,q)+=num[t][q];
            Val(1,q)+=denum[t][q];
        }
    }
}
