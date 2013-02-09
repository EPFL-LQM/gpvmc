#include "SpinSpinCorr.h"
#include "linalg.h"
#include "Stepper.h"
#include "Amplitude.h"
#include "SpinState.h"

SpinSpinCorr::SpinSpinCorr(const Stepper* stepper, FileManager* fm, size_t Nr, size_t* R)
    :MatrixQuantity(stepper,fm,"SpinSpinCorr",Nr,3), m_Nr(Nr)
{
    m_R=new size_t[2*Nr];
    std::memcpy(m_R,R,2*m_Nr*sizeof(size_t));
}

SpinSpinCorr::~SpinSpinCorr() 
{
    delete [] m_R;
}

void SpinSpinCorr::measure()
{
    Quantity::measure();
    const SpinState* st=m_stepper->GetAmp()->GetSpinState();
    for(size_t r=0;r<m_Nr;++r){
        double out(0);
        for(size_t xi=0;xi<st->GetL();++xi){
            for(size_t yi=0;yi<st->GetL();++yi){
                size_t xj=linalg::mod(xi+int(m_R[2*r]),st->GetL());
                size_t yj=linalg::mod(yi+int(m_R[2*r+1]),st->GetL());
                out+=(2*int(st->GetLatOc(xi,yi))-1)*(2*int(st->GetLatOc(xj,yj))-1);
            }
        }
        Val(r,0)=m_R[2*r];
        Val(r,1)=m_R[2*r+1];
        Val(r,2)+=out/std::pow(st->GetL(),2);
    }
}

