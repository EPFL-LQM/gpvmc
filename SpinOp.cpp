#include "SpinOp.h"
#include "Amplitude.h"
#include "SpinState.h"
#include "Jastrow.h"
#include <utility>

using namespace std;

BigComplex SpinOp::Siz(const Amplitude* amp,
                              const size_t& ix, const size_t& iy)
{
    double out;
    const SpinState* st=amp->GetSpinState();
    if(st->GetLatOc(ix,iy)==UP) out=0.5;
    else if(st->GetLatOc(ix,iy)==DOWN) out=-0.5;
    else out=0;
    if(st->GetJas())
        return amp->Amp()*st->GetJas()->Jas()*out;
    else
        return amp->Amp()*out;
}

BigComplex SpinOp::SizSjz(const Amplitude* amp,
                                 const size_t& ix, const size_t& iy,
                                 const size_t& jx, const size_t& jy)
{
    double out;
    const SpinState* st=amp->GetSpinState();
    if((st->GetLatOc(ix,iy)==UP && st->GetLatOc(jx,jy)==UP) || 
       (st->GetLatOc(ix,iy)==DOWN && st->GetLatOc(jx,jy)==DOWN))
        out=0.25;
    else if((st->GetLatOc(ix,iy)==UP && st->GetLatOc(jx,jy)==DOWN) ||
            (st->GetLatOc(ix,iy)==DOWN && st->GetLatOc(jx,jy)==UP))
        out=-0.25;
    else out=0;
    if(st->GetJas())
        return st->GetJas()->Jas()*amp->Amp()*out;
    else
        return amp->Amp()*out;
}


BigComplex SpinOp::SipSjm(const Amplitude* amp,
                          const size_t& ix, const size_t& iy,
                          const size_t& jx, const size_t& jy)
{
    BigComplex out(0);
    const SpinState* st=amp->GetSpinState();
    if(ix==jx && iy==jy){
        if(st->GetLatOc(ix,iy)==UP){
            if(st->GetJas())
                out=st->GetJas()->Jas()*amp->Amp();
            else
                out=amp->Amp();
        } else out=0;
    } else if(st->GetLatOc(ix,iy)==DOWN && st->GetLatOc(jx,jy)==UP){
        size_t idup=st->GetLatUpId(jx,jy), iddo=st->GetLatDoId(ix,iy);
        out=-1.0*amp->VirtRowUpdate(idup,st->GetDo()[iddo],
                                    iddo,st->GetUp()[idup]);
        if(st->GetJas()){
            hop_path_t hopup(1,pair<size_t,size_t>(idup,st->GetDo()[iddo]));
            hop_path_t hopdo(1,pair<size_t,size_t>(iddo,st->GetUp()[idup]));
            out*=st->GetJas()->virtualhop(hopup,hopdo);
        }
    } else {
        out=0;
    }
    return out;
}

BigComplex SpinOp::SimSjp(const Amplitude* amp,
                          const size_t& ix, const size_t& iy,
                          const size_t& jx, const size_t& jy)
{
    const SpinState* st=amp->GetSpinState();
    BigComplex out(0);
    if(ix==jx && iy==jy){
        if(st->GetLatOc(ix,iy)==DOWN){
            if(st->GetJas())
                out=st->GetJas()->Jas()*amp->Amp();
            else
                out=amp->Amp();
        } else out=0;
    } else if(st->GetLatOc(ix,iy)==UP && st->GetLatOc(jx,jy)==DOWN){
        size_t idup=st->GetLatUpId(ix,iy), iddo=st->GetLatDoId(jx,jy);
        out=-1.0*amp->VirtRowUpdate(idup,st->GetDo()[iddo],
                                    iddo,st->GetUp()[idup]);
        if(st->GetJas()){
            hop_path_t hopup(1,pair<size_t,size_t>(idup,st->GetDo()[iddo]));
            hop_path_t hopdo(1,pair<size_t,size_t>(iddo,st->GetUp()[idup]));
            out*=st->GetJas()->virtualhop(hopup,hopdo);
        }
    } else {
        out=0;
    }
    return out;
}
