#include "FullSpaceStepper.h"
#include "SpinState.h"
#include "Jastrow.h"

using namespace std;

void FullSpaceStepper::Reset()
{
    m_prev_up=-1;
    m_prev_down=-1;
    m_khop=-1;
    m_weight=-1;
}

BigDouble FullSpaceStepper::trystep()
{
#ifdef PROFILE
    Timer::tic("FullSpaceStepper::trystep");
#endif
    const WaveFunction* wav=m_amp->GetWaveFunction();
    const SpinState* st=m_amp->GetSpinState();
    size_t idup(st->GetNup()),iddo(st->GetNdo());
    int id1(-1), id2(-1), pick(-1);
    int taus[8]={1,0,-1,0,0,1,0,-1};
    occu_t t1=EMPTY,t2=EMPTY;
    size_t x,y,xp,yp;
    // pick up sites to swap
    x=size_t(st->GetL()*RanGen::uniform());
    y=size_t(st->GetL()*RanGen::uniform());
    t1=st->GetLatOc(x,y);
    if(t1==UP)
        id1=st->GetLatUpId(x,y);
    else if(t1==DOWN)
        id1=st->GetLatDoId(x,y);
    else {
#ifdef EXCEPT
        throw(std::logic_error("FullSpaceStepper::trystep: "
                          "step for empty or double "
                          "occupation "
                          "not implemented"));
#else
        cerr<<"FullSpaceStepper::trystep: "
              "step for empty or double "
              "occupation "
              "not implemented"<<endl;
        abort();
#endif
    }
    pick=int(4.0*RanGen::uniform());
    xp=linalg::mod(x+taus[2*pick],st->GetL());
    yp=linalg::mod(y+taus[2*pick+1],st->GetL());
    t2=st->GetLatOc(xp,yp);
    if(t2==UP)
        id2=st->GetLatUpId(xp,yp);
    else if(t2==DOWN)
        id2=st->GetLatDoId(xp,yp);
    else {
#ifdef EXCEPT
        throw(std::logic_error("FullSpaceStepper::trystep: "
                          "step for empty or double"
                          " occupation "
                          "not implemented"));
#else
        cerr<<"FullSpaceStepper::trystep: "
              "step for empty or double"
              " occupation "
              "not implemented"<<endl;
        abort();
#endif
    }
    if(t1==t2){
        m_prev_up=-1;
        m_prev_down=-1;
        m_khop=-1;
        m_prev=-1;
        return weight();
    } else {
        // done
        // swap sites
        if(t1==UP){
            idup=id1;
            iddo=id2;
        } else if(t1==DOWN){
            idup=id2;
            iddo=id1;
        } else {
#ifdef EXCEPT
            throw(std::logic_error(
                        "FullSpaceStepper::trystep: "
                        "step for empty or double "
                        "occupation not implemented"));
#else
            cerr<<"FullSpaceStepper::trystep: "
                  "step for empty or double "
                  "occupation not implemented"<<endl;
            abort();
#endif
        }
        m_prev_up=idup;
        m_prev_down=iddo;
        vector<hop_path_t> khopup,khopdo,rhopup,rhopdo;
        wav->GetHops(khopup,khopdo);
        rhopup.push_back(
                hop_path_t(1,pair<size_t,size_t>(
                        idup,st->GetDo()[iddo])));
        rhopdo.push_back(
                hop_path_t(1,pair<size_t,size_t>(
                        iddo,st->GetUp()[idup])));
        BigComplex *qs=new BigComplex[khopup.size()];
        BigDouble ampmax(0.0);
#ifdef PROFILE
        Timer::tic("FullSpaceStepper::trystep/VirtUpdate");
#endif
        m_amp->VirtUpdate(rhopup,rhopdo,khopup,khopdo,qs);
#ifdef PROFILE
        Timer::toc("FullSpaceStepper::trystep/VirtUpdate");
#endif
        BigDouble out(0);
        size_t kmax(0);
        for(size_t k=0;k<khopup.size();++k){
                out+=norm(qs[k]);
            if(norm(qs[k])>ampmax){
                ampmax=norm(qs[k]);
                kmax=k;
            }
        }
        m_khop=kmax;
        if(st->GetJas())
            out*=pow(st->GetJas()->virtualhop(rhopup[0],rhopdo[0]),2);
        delete [] qs;
        m_prev=out;
#ifdef PROFILE
        Timer::toc("FullSpaceStepper::trystep");
#endif
        return out;
    }
}

void FullSpaceStepper::step()
{
    SpinState* st=m_amp->GetSpinState();
    WaveFunction *wav=m_amp->GetWaveFunction();
    if(m_prev_up>=0){
        size_t rupid=m_prev_up, rdoid=m_prev_down;
        hop_path_t rhopup(1,
                pair<size_t,size_t>(rupid,st->GetDo()[rdoid]));
        hop_path_t rhopdo(1,
                pair<size_t,size_t>(rdoid,st->GetUp()[rupid]));
        hop_path_t khopup,khopdo;
        wav->GetHop(m_khop,khopup,khopdo);
        st->hop(rhopup,rhopdo);
        wav->hop(m_khop);
        m_amp->Update(rhopup,rhopdo,khopup,khopdo);
    }
    if(m_amp->Amp()==0){
        cerr<<"FullSpaceStepper::step(): Stepped to a zero overlap!"<<endl;
        abort();
    }
    BigDouble ow=m_weight;
#ifndef DEBUG
    if(m_prev>0.0)
        m_weight=m_prev;
    else {
        m_weight=-1;
        weight();
    }
#else
    m_weight=-1;
    weight();
    if(m_prev>=0.0 && abs(double(m_prev/m_weight)-1)>1e-6){
        cerr<<"FullSpaceStepper::step: error: abs(double(m_prev/m_weight))="<<abs(double(m_prev/m_weight))-1<<" m_prev="<<m_prev<<" and m_weight="<<m_weight<<endl;
        abort();
    }
#endif
    if(abs(double(m_weight/ow))<1e-6){
        cerr<<"FullSpaceStepper::step: warning, abs(double(m_weight/ow)) is VERY small! : "<<abs(double(m_weight/ow))<<endl;
    }
    m_prev_up=-1;
    m_prev_down=-1;
    m_khop=-1;
}

BigDouble FullSpaceStepper::weight()
{
    const WaveFunction* wav=m_amp->GetWaveFunction();
    const SpinState* st=m_amp->GetSpinState();
    if(m_weight<0.0){
        BigDouble out(0);
        vector<hop_path_t> khopup,khopdo,rhopup(1),rhopdo(1);
        wav->GetHops(khopup,khopdo);
        vector<BigComplex> qs(khopup.size(),0);
        m_amp->VirtUpdate(rhopup,rhopdo,khopup,khopdo,&qs[0]);
        for(size_t k=0;k<khopup.size();++k){
            out+=norm(qs[k]);
        }
        if(st->GetJas())
            out*=pow(st->GetJas()->Jas(),2);
        m_weight=out;
        if(m_weight==0.0){
            std::cout<<"Warning: FullSpaceStepper::weight: "
                  "state has no overlap, m_prev="
                     <<m_prev<<std::endl;
        }
        return m_weight;
    } else {
        return m_weight;
    }
}

double FullSpaceStepper::transprob()
{
    return 1.0;
}
