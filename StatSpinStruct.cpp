#include <cmath>
#include "StatSpinStruct.h"
#include "Stepper.h"
#include "Jastrow.h"

using namespace std;

StatSpinStruct::StatSpinStruct(const Stepper* stepper,
                               FileManager* fm)
    :MatrixQuantity(stepper,fm,"StatSpinStruct",3,pow(stepper->GetAmp()->GetSpinState()->GetL(),2))
{
    const SpinState* st=stepper->GetAmp()->GetSpinState();
    size_t L=st->GetL();
    m_qs.reserve(L*L);
    for(size_t qx=0; qx<L;++qx){
        for(size_t qy=0;qy<L;++qy){
            m_qs.push_back(qx);
            m_qs.push_back(qy);
        }
    }
    complex<double> I(0,1);
    m_ph=vector<complex<double> >(L*L*m_qs.size()/2);
    for(size_t q=0;q<m_qs.size()/2;++q){
        for(size_t x=0;x<L;++x){
            for(size_t y=0;y<L;++y){
                double phase=2*M_PI*double(m_qs[2*q]*x+m_qs[2*q+1]*y)/L;
                m_ph[(q*L+x)*L+y]=1.0/L*(cos(phase)+I*sin(phase));
            }
        }
    }
}

void StatSpinStruct::measure()
{
    Quantity::measure();
    const SpinState* st=m_stepper->GetAmp()->GetSpinState();
    size_t L=st->GetL();
    // get spin swap list
    vector<hop_path_t> hopup,hopdo;
    for(size_t ix=0;ix<L;++ix){
        for(size_t iy=0;iy<L;++iy){
            for(size_t jx=0;jx<L;++jx){
                for(size_t jy=0;jy<L;++jy){
                    if(st->GetLatOc(ix,iy)==UP && st->GetLatOc(jx,jy)==DOWN){
                        hopup.push_back(
                                hop_path_t(1,hop_t(
                                        st->GetLatUpId(ix,iy),jy*L+jx)));
                        hopdo.push_back(
                                hop_path_t(1,hop_t(
                                        st->GetLatDoId(jx,jy),iy*L+ix)));
                    } else if(st->GetLatOc(jx,jy)==UP && st->GetLatOc(ix,iy)==DOWN){
                        hopup.push_back(
                                hop_path_t(1,hop_t(
                                        st->GetLatUpId(jx,jy),iy*L+ix)));
                        hopdo.push_back(
                                hop_path_t(1,hop_t(
                                        st->GetLatDoId(ix,iy),jy*L+jx)));
                    }
                }
            }
        }
    }
    vector<BigComplex> swamps(hopup.size());
    BigComplex amp=m_stepper->GetAmp()->Amp();
    if(st->GetJas()) amp=amp*st->GetJas()->Jas();
    m_stepper->GetAmp()->VirtUpdate(hopup,hopdo,
                                    vector<hop_path_t>(1),
                                    vector<hop_path_t>(1),
                                    &swamps[0]);
    vector<complex<double> > sqlong(m_qs.size()/2,0);
    vector<BigComplex> sqtranspm(m_qs.size()/2,0);
    vector<BigComplex> sqtransmp(m_qs.size()/2,0);
    for(size_t q=0;q<m_qs.size()/2;++q){
        size_t sw=0;
        for(size_t ix=0;ix<L;++ix){
            for(size_t iy=0;iy<L;++iy){
                for(size_t jx=0;jx<L;++jx){
                    for(size_t jy=0;jy<L;++jy){
                        if(st->GetLatOc(ix,iy)==UP){
                            if(st->GetLatOc(jx,jy)==UP){
                                sqlong[q]+=conj(m_ph[(q*L+jx)*L+jy])*0.25*m_ph[(q*L+ix)*L+iy];
                                if(ix==jx && iy==jy){
                                    sqtranspm[q]+=norm(amp)/double(L*L);
                                }
                            } else if(st->GetLatOc(jx,jy)==DOWN){
                                sqlong[q]+=conj(m_ph[(q*L+jx)*L+jy])*(-0.25)*m_ph[(q*L+ix)*L+iy];
                                if(!st->GetJas()){
                                    sqtransmp[q]-=conj(m_ph[(q*L+jx)*L+jy])*conj(amp)*
                                                  swamps[sw]*m_ph[(q*L+ix)*L+iy];
                                } else {
                                    sqtransmp[q]-=conj(m_ph[(q*L+jx)*L+jy])*conj(amp)*
                                                  swamps[sw]*m_ph[(q*L+ix)*L+iy]*
                                                  st->GetJas()->virtualhop(hopup[sw],hopdo[sw]);
                                }
                                ++sw;
                            }
                        } else if(st->GetLatOc(ix,iy)==DOWN){
                            if(st->GetLatOc(jx,jy)==UP){
                                sqlong[q]+=conj(m_ph[(q*L+jx)*L+jy])*(-0.25)*m_ph[(q*L+ix)*L+iy];
                                if(st->GetJas()){
                                    sqtranspm[q]-=conj(m_ph[(q*L+jx)*L+jy])*conj(amp)*
                                                  swamps[sw]*m_ph[(q*L+ix)*L+iy];
                                } else {
                                    sqtranspm[q]-=conj(m_ph[(q*L+jx)*L+jy])*conj(amp)*
                                                  swamps[sw]*m_ph[(q*L+ix)*L+iy]*
                                                  st->GetJas()->virtualhop(hopup[sw],hopdo[sw]);
                                }
                                ++sw;
                            } else if(st->GetLatOc(jx,jy)==DOWN){
                                sqlong[q]+=conj(m_ph[(q*L+jx)*L+jy])*0.25*m_ph[(q*L+ix)*L+iy];
                                if(ix==jx && iy==jy){
                                    sqtransmp[q]+=norm(amp)/double(L*L);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    BigDouble w=m_stepper->weight();
    for(size_t q=0;q<m_qs.size()/2;++q){
        Val(0,q)+=sqlong[q];
        Val(1,q)+=complex<double>(sqtransmp[q]/w);
        Val(2,q)+=complex<double>(sqtranspm[q]/w);
    }
}
