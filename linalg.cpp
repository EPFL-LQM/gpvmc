#include "linalg.h"
#include "blas_lapack.h"
#include <vector>
#include <algorithm>

using namespace std;

bool linalg::DetInv(const std::complex<double> *A, std::complex<double>* I,
                    size_t M, BigComplex& d)
{
    bool singular(false);
    std::complex<double>* work=new std::complex<double>[M*M];
    memcpy(I,A,M*M*sizeof(std::complex<double>));
    int *ipiv=new int[M], succ, mi(M);
    linalg::zgetrf(mi, mi, I, mi, ipiv ,&succ);
    int perm=1;
    for(size_t i=0;i<M;++i) if(ipiv[i]!=int(i+1)) perm*=-1;
    d=double(perm);
    for(size_t i=0;i<M;++i){
        d*=I[i*M+i];
    }
    singular=norm(d)<1e-6;
    if(!singular)
        linalg::zgetri(mi, I ,mi, ipiv, work, mi, &succ);
    else {
        d=0.0;
#ifdef DEBUG
        std::cout<<"linalg::DetInv: Warning, matrix is singular"<<std::endl;
#endif
    }

    delete [] ipiv;
    delete [] work;

    return singular;
}

bool linalg::Det(const std::complex<double> *A,
                 size_t M, BigComplex& d)
{
    if(!M){
        d=1.0;
        return true;
    } else if(M==1){
        d=A[0];
        return norm(d)<1e-6;
    } else if(M==2){
        d=A[0]*A[3]-A[1]*A[2];
        return norm(d)<1e-6;
    }
    bool singular(false);
    complex<double> *AA=new complex<double>[M*M];
    memcpy(AA,A,M*M*sizeof(complex<double>));
    int *ipiv=new int[M], succ, mi(M);
    linalg::zgetrf(mi, mi, AA, mi, ipiv ,&succ);
    int perm=1;
    for(size_t i=0;i<M;++i) if(ipiv[i]!=int(i+1)) perm*=-1;
    d=double(perm);
    for(size_t i=0;i<M;++i){
        d*=AA[i*M+i];
    }
    singular=norm(d)<1e-6;
    delete [] ipiv;
    delete [] AA;
    return singular;
}

std::string linalg::PrintMat(const std::complex<double>* A,
                             int N, int M, int lda, int precision,bool colwise)
{
    if(!lda) lda=M;
    ostringstream ostr;
    if(colwise){
        for(int i=0;i<N;++i){
            for(int j=0;j<M;++j){
                ostr<<"("<<setw(floor(precision))<<setprecision(precision)
                    <<real(A[j*lda+i]);
                if(imag(A[j*lda+i])>=0)
                    ostr<<"+";
                ostr<<setw(floor(precision))<<setprecision(precision)<<imag(A[j*lda+i])<<"i) ";
            }
            ostr<<endl;
        }
    } else{
        for(int i=0;i<N;++i){
            for(int j=0;j<M;++j){
                ostr<<"("<<setw(floor(precision))<<setprecision(precision)
                    <<real(A[i*lda+j]);
                if(imag(A[i*lda+j])>=0)
                    ostr<<"+";
                ostr<<setw(floor(precision))<<setprecision(precision)<<imag(A[i*lda+j])<<"i) ";
            }
            ostr<<endl;
        }
    }
    return ostr.str();
}

void linalg::DetUpdate(const std::complex<double> *A,
                       const std::complex<double> *Ai,
                       size_t N,
                       const std::complex<double> *C,
                       size_t ldc,
                       const size_t * ci,
                       const size_t * rc,
                       const size_t nc,
                       const std::complex<double> *R,
                       size_t ldr,
                       const size_t * ri,
                       const size_t * rr,
                       const size_t nr,
                       BigComplex *dets)
{
#ifdef PROFILE
    Timer::tic("DetUpdate");
#endif
    // vectors of indices for heterogeneous
    // rows/columns ranks
    vector<size_t> rid(nr,0),cid(nc,0);
    for(size_t r=1;r<nr;++r) rid[r]=rid[r-1]+rr[r-1];
    for(size_t c=1;c<nc;++c) cid[c]=cid[c-1]+rc[c-1];
    // total number of rows/columns
    size_t Nc=cid.back()+rc[nc-1];
    size_t Nr=rid.back()+rr[nr-1];
    // max rank for rows/columns
    size_t rmax=0;
    size_t cmax=0;
    for(size_t r=0;r<nr;++r) if(rr[r]>rmax) rmax=rr[r];
    for(size_t c=0;c<nc;++c) if(rc[c]>cmax) cmax=rc[c];
    if(nc<=1){
        complex<double> *EcAiC=new complex<double>[rc[0]*rc[0]];
        complex<double> *AiC=new complex<double>[N*rc[0]];
        complex<double> *RAiC=new complex<double>[Nr*rc[0]];
        complex<double> one(1,0), zer(0,0), mone(-1,0);
        if(nc && rc[0]){
            //calculate Ai*C
            if(rc[0]==1)
                linalg::zgemv('N',N,N,&one,Ai,N,C,1,&zer,AiC,1);
            else
                linalg::zgemm('N','N',N,rc[0],N,
                              &one,Ai,N,C,ldc,&zer,AiC,N);
            //cout<<"AiC="<<endl<<PrintMat(AiC,N,rc,N)<<endl;
            //copy Ec*Ai*C
            for(size_t i=0;i<rc[0];++i)
                for(size_t j=0;j<rc[0];++j)
                    EcAiC[j*rc[0]+i]=AiC[ci[i]+j*N];
            //cout<<"EcAiC="<<endl<<PrintMat(EcAiC,rc,rc,rc)<<endl;
            //calculate R*Ai*C
            if(Nr){
                if(rc[0]==1)
                    linalg::zgemv('N',Nr,N,&one,R,ldr,
                                  AiC,1,&zer,RAiC,1);
                else
                    linalg::zgemm('N','N',Nr,rc[0],N,
                                  &one,R,ldr,AiC,N,&zer,RAiC,Nr);
            }
        }
        complex<double> *B=new complex<double>[(cmax+rmax)*(cmax+rmax)];
        complex<double> *ErTAEcT=new complex<double>[rmax*cmax];
        complex<double> *ErTC=new complex<double>[rmax*cmax];
        for(size_t r=0;r<nr;++r){
            // calculate R*Ai*Er and fill it in B22
            for(size_t i=0;i<rr[r];++i){
                for(size_t j=0;j<rr[r];++j){
                    //cout<<"R["<<r*rr+i<<"]="<<PrintMat(&R[r*rr+i],1,N,rr*nr)<<endl;
                    //cout<<"Ai["<<ri[r*rr+j]<<"]="<<PrintMat(&Ai[ri[r*rr+j]*N],N,1,N)<<endl;
                    linalg::zdotu_sub(N,&R[rid[r]+i],ldr,
                                      &Ai[ri[rid[r]+j]*N],1,
                                      &B[rc[0]+i+(rc[0]+j)*(rr[r]+rc[0])]);
                }
            }
            // fill B11 with simple elements
            for(size_t i=0;i<rc[0];++i)
                for(size_t j=0;j<rc[0];++j)
                    B[j*(rc[0]+rr[r])+i]=EcAiC[j*rc[0]+i];
            //std::cout<<"B22="<<std::endl<<linalg::PrintMat(&B[rc+rc*(rr+rc)],rr,rr,rr+rc)<<std::endl;
            if(nc && rc[0] && rr[r]){
                // copy Ec*Ai*Er in top right block in B
                for(size_t i=0;i<rc[0];++i)
                    for(size_t j=0;j<rr[r];++j)
                        B[(rc[0]+j)*(rc[0]+rr[r])+i]=Ai[ci[i]+ri[rid[r]+j]*N];
                //std::cout<<"B12="<<std::endl<<PrintMat(&B[rc*(rr+rc)],rc,rr,rc+rr)<<std::endl;
                // copy ErT*A*EcT
                // copy ErT*C
                // fill simple elements in B21
                for(size_t i=0;i<rr[r];++i){
                    for(size_t j=0;j<rc[0];++j){
                        ErTAEcT[j*rr[r]+i]=A[ri[rid[r]+i]+ci[j]*N];
                        ErTC[j*rr[r]+i]=C[ri[rid[r]+i]+j*N];
                        B[rc[0]+i+j*(rr[r]+rc[0])]=RAiC[rid[r]+i+j*Nr]
                                                 -R[rid[r]+i+ci[j]*ldr];
                    }
                }
                //std::cout<<"ErTAEcT="<<std::endl<<PrintMat(ErTAEcT,rr,rc,rr)<<endl;
                //cout<<"ErTC="<<endl<<PrintMat(ErTC,rr,rc,rr)<<endl;
                //cout<<"B21="<<endl<<PrintMat(&B[rc],rr,rc,rr+rc)<<endl;
                // fill B21 with multiplicative elements
                if(rr[r]==1){
                    if(rc[0]==1){
                        B[rc[0]]+=B[rc[0]*(rr[r]+rc[0])+rc[0]]*(ErTAEcT[0]+ErTC[0]);
                    } else {
                        for(size_t i=0;i<rc[0];++i)
                            B[rc[0]+i*(rr[r]+rc[0])]+=B[rc[0]*(rr[r]+rc[0])+rc[0]]*(ErTAEcT[i]+ErTC[i]);
                    }
                } else if(rc[0]==1) {
                    linalg::zgemv('N',rr[r],rr[r],&one,&B[rc[0]*(rr[r]+rc[0])+rc[0]],rc[0]+rr[r],
                                  ErTAEcT,1,&one,&B[rc[0]],1);
                    //cout<<"B21p="<<endl<<PrintMat(&B[rc],rr,rc,rr+rc)<<endl;
                    linalg::zgemv('N',rr[r],rr[r],&one,&B[rc[0]*(rr[r]+rc[0])+rc[0]],rc[0]+rr[r],
                                  ErTC,1,&one,&B[rc[0]],1);
                } else {
                    linalg::zgemm('N','N',rr[r],rc[0],rr[r],
                                  &one,&B[rc[0]*(rr[r]+rc[0])+rc[0]],rr[r]+rc[0],ErTAEcT,rr[r],
                                  &one,&B[rc[0]],rc[0]+rr[r]);
                    linalg::zgemm('N','N',rr[r],rc[0],rr[r],
                                  &mone,&B[rc[0]*(rr[r]+rc[0])+rc[0]],rr[r]+rc[0],ErTC,rr[r],
                                  &one,&B[rc[0]],rc[0]+rr[r]);
                }
                //cout<<"B21="<<endl<<PrintMat(&B[rc],rr,rc,rr+rc)<<endl;
                //cout<<"B11="<<endl<<PrintMat(B,rc,rc,rr+rc)<<endl;
                // add multiplicative elements to B11
                if(rr[r]==1 && rc[0]==1){
                    B[0]+=B[rc[0]+rr[r]]*(ErTC[0]+ErTAEcT[0]);
                } else if(rr[r]==1){
                    //cout<<"oper is"<<endl<<PrintMat(&B[rc*(rr+rc)],rc,1,rc)
                        //<<endl<<"times"<<endl
                        //<<PrintMat(ErTC,1,rc,1)<<endl;
                    linalg::zgeru(rc[0],rc[0],&one,&B[rc[0]*(rc[0]+rr[r])],1,ErTC,1,
                                  B,rc[0]+rr[r]);
                    //cout<<"B11="<<endl<<PrintMat(B,rc,rc,rr+rc)<<endl;
                    linalg::zgeru(rc[0],rc[0],&one,&B[rc[0]*(rc[0]+rr[r])],1,ErTAEcT,1,
                                  B,rc[0]+rr[r]);
                } else if(rc[0]==1){
                    complex<double> mem;
                    linalg::zdotu_sub(rr[r],&B[rc[0]*(rc[0]+rr[r])],rc[0]+rr[r],ErTC,1,&mem);
                    B[0]+=mem;
                    linalg::zdotu_sub(rr[r],&B[rc[0]*(rc[0]+rr[r])],rc[0]+rr[r],ErTAEcT,1,&mem);
                    B[0]+=mem;
                } else {
                    linalg::zgemm('N','N',rc[0],rc[0],rr[r],
                                  &mone,&B[rc[0]*(rc[0]+rr[r])],rc[0]+rr[r],ErTC,rr[r],
                                  &one,B,rc[0]+rr[r]);
                    linalg::zgemm('N', 'N',rc[0],rc[0],rr[r],
                                  &one,&B[rc[0]*(rc[0]+rr[r])],rc[0]+rr[r],ErTAEcT,rr[r],
                                  &one,B,rc[0]+rr[r]);
                }
            }
            //cout<<"B="<<endl<<PrintMat(B,rc+rr,rc+rr,rc+rr)<<endl;
            linalg::Det(B,rc[0]+rr[r],dets[r]);
            //cout<<"det="<<dets[r]<<endl;
        }
        delete [] ErTAEcT;
        delete [] ErTC;
        delete [] B;
        delete [] RAiC;
        delete [] AiC;
        delete [] EcAiC;
    } else if(nr<=1){
        complex<double> *RAiEr=new complex<double>[rr[0]*rr[0]];
        complex<double> *RAi=new complex<double>[rr[0]*N];
        complex<double> *RAiC=new complex<double>[rr[0]*Nc];
        complex<double> one(1,0), mone(-1,0), zer(0,0);
        if(nr && rr[0]){
            // calculate R*Ai
            if(rr[0]==1)
                linalg::zgemv('T',N,N,&one,Ai,N,R,ldr,&zer,RAi,1);
            else
                linalg::zgemm('N','N',rr[0],N,N,
                              &one,R,ldr,Ai,N,&zer,RAi,rr[0]);
            // copy RAiEr
            for(size_t i=0;i<rr[0];++i)
                for(size_t j=0;j<rr[0];++j)
                    RAiEr[j*rr[0]+i]=RAi[i+ri[j]*rr[0]];
            //cout<<"RAiEr="<<endl<<PrintMat(RAiEr,rr,rr,rr)<<endl;
            // calculate RAiC
            if(Nc){
                if(rr[0]==1)
                    linalg::zgemv('T',N,Nc,&one,C,ldc,
                                  RAi,1,&zer,RAiC,1);
                else
                    linalg::zgemm('N','N',rr[0],Nc,N,
                                  &one,RAi,rr[0],C,ldc,&zer,RAiC,rr[0]);
            }
            //cout<<"RAiC="<<endl<<PrintMat(RAiC,rr,nc*rc,rr)<<endl;
        }
        complex<double> *B=new complex<double>[(cmax+rmax)*(cmax+rmax)];
        complex<double> *ErTAEcT=new complex<double>[rmax*cmax];
        complex<double> *REcT=new complex<double>[rmax*cmax];
        for(size_t c=0;c<nc;++c){
            // calculate EcAiC and fill B11
            for(size_t i=0;i<rc[c];++i)
                for(size_t j=0;j<rc[c];++j)
                    linalg::zdotu_sub(N,&Ai[ci[cid[c]+i]],N,
                                      &C[(cid[c]+j)*ldc],1,
                                      &B[j*(rc[c]+rr[0])+i]);
            // copy R*Ai*Er in B22
            for(size_t i=0;i<rr[0];++i)
                for(size_t j=0;j<rr[0];++j)
                    B[(rc[c]+j)*(rc[c]+rr[0])+rc[c]+i]=RAiEr[j*rr[0]+i];
            //cout<<"B11="<<endl<<PrintMat(B,rc,rc,rc+rr)<<endl;
            if(rr[0] && rc[c]){
                // copy Ec*Ai*Er in B12
                for(size_t i=0;i<rc[c];++i)
                    for(size_t j=0;j<rr[0];++j)
                        B[i+(rc[c]+j)*(rc[c]+rr[0])]=Ai[ci[cid[c]+i]+ri[j]*N];
                //cout<<"B12="<<endl<<PrintMat(&B[rc*(rc+rr)],rc,rr,rc+rr)<<endl;
                // copy ErT*A*EcT
                // copy REcT
                // fill simple elements in B21
                for(size_t i=0;i<rr[0];++i){
                    for(size_t j=0;j<rc[c];++j){
                        ErTAEcT[j*rr[0]+i]=A[ri[i]+ci[cid[c]+j]*N];
                        REcT[j*rr[0]+i]=R[i+ci[cid[c]+j]*nr*rr[0]];
                        B[rc[c]+i+j*(rc[c]+rr[0])]=RAiC[i+(j+cid[c])*rr[0]]
                                            -C[ri[i]+(cid[c]+j)*ldc];
                    }
                }
                //cout<<"RAiC[n]="<<endl<<PrintMat(&RAiC[c*rc*rr],rr,rc,rr)<<endl;
                //cout<<"ErTAEcT="<<endl<<PrintMat(ErTAEcT,rr,rc,rr)<<endl;
                //cout<<"REcT="<<endl<<PrintMat(REcT,rr,rc,rr)<<endl;
                //cout<<"B21="<<endl<<PrintMat(&B[rc],rr,rc,rc+rr)<<endl;
                // fill B21 with multiplicative elements
                if(rc[c]==1){
                    if(rr[0]==1){
                        B[rc[c]]+=(ErTAEcT[0]+REcT[0])*B[0];
                    } else {
                        for(size_t i=0;i<rr[0];++i)
                            B[(rc[c]+i)]+=(ErTAEcT[i]+REcT[i])*B[0];
                    }
                } else if(rr[0]==1){
                    linalg::zgemv('T',rc[c],rc[c],&one,B,rc[c]+rr[0],
                                  ErTAEcT,1,&one,&B[rc[c]],rc[c]+rr[0]);
                    linalg::zgemv('T',rc[c],rc[c],&one,B,rc[c]+rr[0],
                                  REcT,1,&one,&B[rc[c]],rc[c]+rr[0]);
                } else {
                    linalg::zgemm('N','N',rr[0],rc[c],rc[c],
                                  &one,ErTAEcT,rr[0],B,rc[c]+rr[0],
                                  &one,&B[rc[c]],rc[c]+rr[0]);
                    linalg::zgemm('N','N',rr[0],rc[c],rc[c],
                                  &mone,REcT,rr[0],B,rc[c]+rr[0],
                                  &one,&B[rc[c]],rc[c]+rr[0]);
                }
                //cout<<"B21="<<endl<<PrintMat(&B[rc],rr,rc,rr+rc)<<endl;
                // add multiplicative elements to B22
                if(rc[c]==1 && rr[0]==1){
                    B[rc[c]*(rc[c]+rr[0])+rc[c]]+=(REcT[0]+ErTAEcT[0])*B[rc[c]*(rr[0]+rc[c])];
                } else if(rc[c]==1){
                    linalg::zgeru(rr[0],rr[0],&one,REcT,1,&B[rc[c]*(rr[0]+rc[c])],rr[0]+rc[c],
                                  &B[rc[c]*(rc[c]+rr[0])+rc[c]],rc[c]+rr[0]);
                    linalg::zgeru(rr[0],rr[0],&one,ErTAEcT,1,&B[rc[c]*(rr[0]+rc[c])],rr[0]+rc[c],
                                  &B[rc[c]*(rc[c]+rr[0])+rc[c]],rc[c]+rr[0]);
                } else if(rr[0]==1){
                    complex<double> mem;
                    linalg::zdotu_sub(rc[c],REcT,1,&B[rc[c]*(rr[0]+rc[c])],1,&mem);
                    B[rc[c]*(rc[c]+rr[0])+rc[c]]+=mem;
                    linalg::zdotu_sub(rc[c],ErTAEcT,1,&B[rc[c]*(rr[0]+rc[c])],1,&mem);
                    B[rc[c]*(rc[c]+rr[0])+rc[c]]+=mem;
                } else {
                    linalg::zgemm('N','N',rr[0],rr[0],rc[c],
                                  &mone,REcT,rr[0],&B[rc[c]*(rr[0]+rc[c])],rc[c]+rr[0],
                                  &one,&B[rc[c]*(rc[c]+rr[0])+rc[c]],rc[c]+rr[0]);
                    linalg::zgemm('N','N',rr[0],rr[0],rc[c],
                                  &one,ErTAEcT,rr[0],&B[rc[c]*(rr[0]+rc[c])],rc[c]+rr[0],
                                  &one,&B[rc[c]*(rc[c]+rr[0])+rc[c]],rc[c]+rr[0]);
                }
            }
            //cout<<"B="<<endl<<PrintMat(B,rr+rc,rr+rc,rr+rc)<<endl;
            linalg::Det(B,rc[c]+rr[0],dets[c]);
            //cout<<"det="<<dets[c]<<endl;
        }
        delete [] ErTAEcT;
        delete [] REcT;
        delete [] B;
        delete [] RAiEr;
        delete [] RAi;
        delete [] RAiC;
    } else {
#ifdef EXCEPT
        throw(std::logic_error("linalg::DetUpdate: "
                               "it is required that "
                               "min(nc,nr)=1"));
#else
        abort();
#endif
    }
#ifdef PROFILE
    Timer::toc("DetUpdate");
#endif
}

void linalg::InvUpdate(const std::complex<double> *A,
                       std::complex<double> *Ai,
                       size_t N,
                       const std::complex<double> *C,
                       const size_t *ci,
                       size_t rc,
                       const std::complex<double> *R,
                       const size_t *ri,
                       size_t rr, BigComplex& det)
{
    if(!rc && !rr){
        det=1;
        return;
    }
    complex<double> *buf=new complex<double>[3*rc*rc+3*rr*rr+
                                             8*rc*rr+2*N*(rc+rr)+
                                             max(rc,rr)*N];
    complex<double> *EcAiC=  &buf[0];//[rc*rc];
    complex<double> *RAiEr=  &buf[rc*rc];//[rr*rr];
    complex<double> *B=      &buf[rc*rc+rr*rr];//[(rc+rr)*(rc+rr)];
    complex<double> *Bi=     &buf[2*rc*rc+2*rr*rr+2*rc*rr];//[(rc+rr)*(rc+rr)];
    complex<double> *RAiC=   &buf[3*rc*rc+3*rr*rr+4*rc*rr];//[rr*rc];
    complex<double> *ErTAEcT=&buf[3*rc*rc+3*rr*rr+5*rc*rr];//[rr*rc];
    complex<double> *ErTC=   &buf[3*rc*rc+3*rr*rr+6*rc*rr];//[rr*rc];
    complex<double> *REcT=   &buf[3*rc*rc+3*rr*rr+7*rc*rr];//[rr*rc];
    complex<double> *AiC=    &buf[3*rc*rc+3*rr*rr+8*rc*rr];//[N*rc];
    complex<double> *EcAi=   &buf[3*rc*rc+3*rr*rr+8*rc*rr+N*rc];//[rc*N];
    complex<double> *RAi=    &buf[3*rc*rc+3*rr*rr+8*rc*rr+2*N*rc];//[rr*N];
    complex<double> *AiEr=   &buf[3*rc*rc+3*rr*rr+8*rc*rr+N*(2*rc+rr)];//[N*rr]
    complex<double> *mem=    &buf[3*rc*rc+3*rr*rr+8*rc*rr+2*N*(rc+rr)];//[max(rr,rc)*N]
    complex<double> one(1,0), zer(0,0), mone(-1,0);
    //calculate Ai*C
    if(rc) linalg::zgemm('N','N',N,rc,N,
                         &one,Ai,N,C,N,&zer,AiC,N);
    // calculate R*Ai
    if(rr) linalg::zgemm('N','N',rr,N,N,
                         &one,R,rr,Ai,N,
                         &zer,RAi,rr);
    //copy Ec*Ai*C
    for(size_t i=0;i<rc;++i)
        for(size_t j=0;j<rc;++j)
            EcAiC[j*rc+i]=AiC[ci[i]+j*N];
    // copy R*Ai*Er
    for(size_t i=0;i<rr;++i)
        for(size_t j=0;j<rr;++j)
            RAiEr[i+j*rr]=RAi[i+ri[j]*rr];
    // calculate RAiC
    if(rc*rr) linalg::zgemm('N','N',rr,rc,N,
                            &one,RAi,rr,C,N,
                            &zer,RAiC,rr);
    // copy ErTAEcT
    // copy REcT
    // copy ErTC
    for(size_t i=0;i<rr;++i){
        for(size_t j=0;j<rc;++j){
            ErTAEcT[i+j*rr]=A[ri[i]+ci[j]*N];
            ErTC[i+j*rr]=C[ri[i]+j*N];
            REcT[i+j*rr]=R[i+ci[j]*rr];
        }
    }
    // copy Ec*Ai*Er in B12
    for(size_t i=0;i<rc;++i)
        for(size_t j=0;j<rr;++j)
            B[i+(rc+j)*(rc+rr)]=Ai[ci[i]+ri[j]*N];
    //cout<<"B12="<<endl<<PrintMat(&B[rc*(rc+rr)],rc,rr,rc+rr)<<endl;
    // copy R*Ai*Er in B22
    for(size_t i=0;i<rr;++i)
        for(size_t j=0;j<rr;++j)
            B[(rc+i)+(rc+j)*(rc+rr)]=RAiEr[i+j*rr];
    // fill B11
    for(size_t i=0;i<rc;++i)
        for(size_t j=0;j<rc;++j)
            B[i+j*(rc+rr)]=EcAiC[i+j*rc];
    //cout<<"B11="<<endl<<PrintMat(B,rc,rc,rr+rc)<<endl;
    if(rr*rc){
        linalg::zgemm('N','N',rc,rc,rr,
                      &mone,&B[rc*(rc+rr)],rc+rr,ErTC,rr,
                      &one,B,rc+rr);
        linalg::zgemm('N','N',rc,rc,rr,
                      &one,&B[rc*(rc+rr)],rc+rr,ErTAEcT,rr,
                      &one,B,rc+rr);
    }
    // fill B21
    for(size_t i=0;i<rr;++i){
        for(size_t j=0;j<rc;++j){
            B[(rc+i)+j*(rc+rr)]=RAiC[i+j*rr]
                                -R[i+ci[j]*rr];
        }
    }
    if(rr*rc){
        linalg::zgemm('N','N',rr,rc,rr,
                      &one,&B[rc*(rr+rc)+rc],rr+rc,ErTAEcT,rr,
                      &one,&B[rc],rc+rr);
        linalg::zgemm('N','N',rr,rc,rr,
                      &mone,&B[rc*(rr+rc)+rc],rr+rc,ErTC,rr,
                      &one,&B[rc],rc+rr);
    }
    //cout<<"B="<<PrintMat(B,rc+rr,rc+rr,rc+rr)<<endl;
    if(!linalg::DetInv(B,Bi,(rr+rc),det)){
        //cout<<"det="<<d<<endl;
        // copy AiEr ans EcAi
        for(size_t i=0;i<N;++i){
            for(size_t j=0;j<rr;++j)
                AiEr[i+j*N]=Ai[i+ri[j]*N];
            for(size_t j=0;j<rc;++j)
                EcAi[j+i*rc]=Ai[ci[j]+i*N];
        }
        // update AiC so it matches Ai*Cp=Ai*(1-Er*ErT)*(C-A*EcT)
        if(rr*rc){
            linalg::zgemm('N','N',N,rc,rr,
                          &mone,AiEr,N,ErTC,rr,
                          &one,AiC,N);
            linalg::zgemm('N','N',N,rc,rr,
                          &one,AiEr,N,ErTAEcT,rr,
                          &one,AiC,N);
        }
        for(size_t c=0;c<rc;++c)
            AiC[ci[c]+c*N]-=1;
        // update RAi sot it matches Rp*Ai=(R-ErT*A)*Ai
        for(size_t r=0;r<rr;++r)
            RAi[r+ri[r]*rr]-=1;
        // Update inverse matrix
        if(rc){
            // use mem to store B11*Ec*Ai
            linalg::zgemm('N','N',rc,N,rc,
                          &one,Bi,rr+rc,EcAi,rc,
                          &zer,mem,rc);
            // add -Ai*Cp*B11*Ec*Ai in Ai
            linalg::zgemm('N','N',N,N,rc,
                          &mone,AiC,N,mem,rc,
                          &one,Ai,N);
        }
        if(rc*rr){
            // use mem to store B12*Rp*Ai
            linalg::zgemm('N','N',rc,N,rr,
                          &one,&Bi[rc*(rc+rr)],rr+rc,RAi,rr,
                          &zer,mem,rc);
            // add -Ai*Cp*B12*Rp*Ai to in Ai
            linalg::zgemm('N','N',N,N,rc,
                          &mone,AiC,N,mem,rc,
                          &one,Ai,N);
            // use mem to store B21*Ec*Ai
            linalg::zgemm('N','N',rr,N,rc,
                          &one,&Bi[rc],rr+rc,EcAi,rc,
                          &zer,mem,rr);
            // add -Ai*Er*B21*Ec*Ai to Ai
            linalg::zgemm('N','N',N,N,rr,
                          &mone,AiEr,N,mem,rr,
                          &one,Ai,N);
        }
        if(rr){
            // use mem to store B22*Rp*Ai
            linalg::zgemm('N','N',rr,N,rr,
                          &one,&Bi[rc*(rr+rc)+rc],rc+rr,RAi,rr,
                          &zer,mem,rr);
            // add -Ai*Er*B22*Rp*Ai to Ai
            linalg::zgemm('N','N',N,N,rr,
                          &mone,AiEr,N,mem,rr,
                          &one,Ai,N);
        }
    } else {
#ifdef EXCEPT
        throw(std::logic_error("linalg::InvUpdate: Updated matrix nearly "
                               "singular."));
#else
        abort();
#endif
    }
    delete [] buf;
}
