#include "HeisEn.h"
#include "SpinOp.h"
#include "Stepper.h"
#include "Amplitude.h"
#include "SpinState.h"
#include "linalg.h"
#include "Jastrow.h"

void HeisEn::measure()
{
    const Amplitude* amp=m_stepper->GetAmp();
    Quantity::measure();
    const SpinState* st=amp->GetSpinState();
    BigComplex out(0);
    int taus[16]={1,0,-1,0,0,1,0,-1,1,1,-1,-1,1,-1,-1,1};
    double Js[8]={1,1,1,1,m_jr,m_jr,m_jr,m_jr};
    size_t nm=4;
    if(m_jr!=0) nm=8;
    for(size_t x=0;x<st->GetL();++x){
        for(size_t y=0;y<st->GetL();++y){
            for(size_t t=0;t<nm;++t){
                size_t jx=linalg::mod(int(x)+taus[2*t],st->GetL());
                size_t jy=linalg::mod(int(y)+taus[2*t+1],st->GetL());
                out+=Js[t]*SpinOp::SizSjz(amp,x,y,jx,jy)+0.5*(SpinOp::SipSjm(amp,x,y,jx,jy)+SpinOp::SimSjp(amp,x,y,jx,jy));
            }
        }
    }
    out*=conj(amp->Amp());
    if(st->GetJas())
        out*=st->GetJas()->Jas();
    BigDouble w=m_stepper->weight();
    m_vals+=std::complex<double>(0.5*out/double(st->GetL()*st->GetL())/w);
}
