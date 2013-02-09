#include "SpinState.h"
#include "Jastrow.h"
#include "linalg.h"

using namespace std;

SpinState::SpinState(size_t L, size_t Nup, size_t Ndown, bool neel, bool doccu)
    :m_L(L), m_Nup(Nup), m_Ndo(Ndown), m_jas(0), m_own_amp(false), m_sign(1)
{
#ifdef DEBUG
    if((!doccu && m_Nup+m_Ndo>m_L*m_L) ||
        (doccu && m_Nup+m_Ndo>2*m_L*m_L))
#ifdef EXCEPT
        throw(std::out_of_range("SpinState::SpinState"
                           "(size_t, size_t, size_t,bool): "
                           "Too many fermions to fill the lattice."));
#else
        abort();
#endif
#endif
    m_up=new size_t[m_Nup];
    m_do=new size_t[m_Ndo];
    m_latupid=new size_t[m_L*m_L];
    m_latdoid=new size_t[m_L*m_L];
    m_latoc=new occu_t[m_L*m_L];

    Init(neel,doccu);    
}

void SpinState::Init(bool neel, bool doccu)
{
    for(size_t x=0;x<m_L;++x){
        for(size_t y=0;y<m_L;++y){
            m_latoc[y*m_L+x]=EMPTY;
            m_latupid[y*m_L+x]=m_Nup;
            m_latdoid[y*m_L+x]=m_Ndo;
        }
    }
    if(neel){
        if(m_L%2!=0)
#ifdef EXCEPT
            throw(std::logic_error("SpinState::SpinState"
                              "(size_t, size_t, size_t,bool): "
                              "Neel state require an even "
                              "number of sites."));
#else
            abort();
#endif
        size_t nup=0,ndown=0;
        std::vector<size_t> flip(0,m_L*m_L);
        int tf=int(m_Nup)-int(m_Ndo);
        while(flip.size()!=size_t(abs(tf)/2)){
            size_t pick=floor(RanGen::uniform()*m_L*m_L/2);
            while(find(flip.begin(),flip.end(),pick)!=flip.end())
                pick=floor(RanGen::uniform()*m_L*m_L/2);
            flip.push_back(pick);
        }
        for(size_t x=0;x<m_L;++x){
            for(size_t y=0;y<m_L;++y){
                if((x+y)%2!=0){
                    if(tf<0 && find(flip.begin(),
                                    flip.end(),(x*m_L+y)/2)!=flip.end()){
                        m_latoc[y*m_L+x]=DOWN;
                        m_latdoid[y*m_L+x]=ndown;
                        m_do[ndown]=y*m_L+x;
                        m_sign*=cdagsign(m_do[ndown],false);
                        ++ndown;
                    } else {
                        m_latoc[y*m_L+x]=UP;
                        m_latupid[y*m_L+x]=nup;
                        m_up[nup]=y*m_L+x;
                        m_sign*=cdagsign(m_up[nup],true);
                        ++nup;
                    }
                } else {
                    if(tf>0 && find(flip.begin(),
                                    flip.end(),(x*m_L+y)/2)!=flip.end()){
                        m_latoc[y*m_L+x]=UP;
                        m_latupid[y*m_L+x]=nup;
                        m_up[nup]=y*m_L+x;
                        m_sign*=cdagsign(m_up[nup],true);
                        ++nup;
                    } else {
                        m_latoc[y*m_L+x]=DOWN;
                        m_latdoid[y*m_L+x]=ndown;
                        m_do[ndown]=y*m_L+x;
                        m_sign*=cdagsign(m_do[ndown],false);
                        ++ndown;
                    }
                }
            }
        }
    } else {
        // place randomly spins up.
        size_t upleft=m_Nup;
        while(upleft!=0){
            size_t x=floor(RanGen::uniform()*m_L);
            size_t y=floor(RanGen::uniform()*m_L);
            if(m_latoc[y*m_L+x]==EMPTY){
                m_latoc[y*m_L+x]=UP;
                m_latupid[y*m_L+x]=size_t(int(m_Nup)-int(upleft));
                m_up[(int(m_Nup)-int(upleft))]=y*m_L+x;
                m_sign*=cdagsign(y*m_L+x,true);
                upleft--;
            }
        }

        // place randomly spins down.
        size_t downleft=m_Ndo;
        while(downleft!=0){
            size_t x=floor(RanGen::uniform()*m_L);
            size_t y=floor(RanGen::uniform()*m_L);
            if(m_latoc[y*m_L+x]==EMPTY){
                m_latoc[y*m_L+x]=DOWN;
                m_latdoid[y*m_L+x]=size_t(int(m_Ndo)-int(downleft));
                m_do[(int(m_Ndo)-int(downleft))]=y*m_L+x;
                downleft--;
                m_sign*=cdagsign(y*m_L+x,false);
            } else if(doccu && m_latoc[y*m_L+x]==UP){
                m_latoc[y*m_L+x]=DO;
                m_latdoid[y*m_L+x]=size_t(int(m_Ndo)-int(downleft));
                m_do[(int(m_Ndo)-int(downleft))]=y*m_L+x;
                downleft--;
            }
        }
    }
}

void SpinState::Init(size_t * state)
{
    for(size_t x=0;x<m_L;++x){
        for(size_t y=0;y<m_L;++y){
            m_latoc[y*m_L+x]=EMPTY;
            m_latupid[y*m_L+x]=m_Nup;
            m_latdoid[y*m_L+x]=m_Ndo;
        }
    }
    size_t ui(0),di(0);
    for(size_t x=0; x<m_L; ++x){
        for(size_t y=0;y<m_L;++y){
            if(state[x*m_L+y]==1){
                m_latoc[x*m_L+y]=UP;
                m_latupid[x*m_L+y]=ui;
                m_up[ui]=y*m_L+x;
                ++ui;
            } else {
                m_latoc[x*m_L+y]=DOWN;
                m_latdoid[x*m_L+y]=di;
                m_do[di]=y*m_L+x;
                ++di;
            }
        }
    }
    //m_nn=count_nn();
}

SpinState::~SpinState()
{
    delete [] m_up;
    delete [] m_do;
    delete [] m_latupid;
    delete [] m_latdoid;
    delete [] m_latoc;
}

void SpinState::hop(const hop_path_t& hopup,
                    const hop_path_t& hopdo)
{
    //swap two spins and update the corresponding amplitude
    if(m_jas) m_jas->hop(hopup,hopdo);
    m_sign*=hop_sign(hopup,hopdo);
    for(size_t h=0;h<hopup.size();++h){
        m_latupid[m_up[hopup[h].first]]=m_Nup;
        m_latupid[hopup[h].second]=hopup[h].first;
        if(m_latoc[m_up[hopup[h].first]]==UP)
            m_latoc[m_up[hopup[h].first]]=EMPTY;
        else
            m_latoc[m_up[hopup[h].first]]=DOWN;
        if(m_latoc[hopup[h].second]==DOWN)
            m_latoc[hopup[h].second]=DO;
        else
            m_latoc[hopup[h].second]=UP;
        m_up[hopup[h].first]=hopup[h].second;
    }
    for(size_t h=0;h<hopdo.size();++h){
        m_latdoid[m_do[hopdo[h].first]]=m_Ndo;
        m_latdoid[hopdo[h].second]=hopdo[h].first;
        if(m_latoc[m_do[hopdo[h].first]]==DOWN)
            m_latoc[m_do[hopdo[h].first]]=EMPTY;
        else
            m_latoc[m_do[hopdo[h].first]]=UP;
        if(m_latoc[hopdo[h].second]==UP)
            m_latoc[hopdo[h].second]=DO;
        else
            m_latoc[hopdo[h].second]=DOWN;
        m_do[hopdo[h].first]=hopdo[h].second;
    }
}

double SpinState::Jas() const
{
    if(m_jas)
        return m_jas->Jas();
    else
        return 1.0;
}

int SpinState::hop_sign(const hop_path_t& hopup,
                        const hop_path_t& hopdo) const
{
    bool odd=false;
    size_t mn,mx;
    if(hopup.size()){
        mn=min(hopup[0].second,m_up[hopup[0].first])+1;
        mx=max(hopup[0].second,m_up[hopup[0].first]);
        for(size_t s=mn;s<mx;++s)
            odd= (odd != (m_latupid[s]!=m_Nup));
        if(hopup.size()>1){
            size_t* sup=new size_t[m_L*m_L];
            memcpy(sup,m_latupid,m_L*m_L*sizeof(size_t));
            sup[m_up[hopup[0].first]]=m_Nup;
            sup[hopup[0].second]=hopup[0].first;
            for(size_t h=1;h<hopup.size();++h){
                mn=min(hopup[0].second,m_up[hopup[0].first])+1;
                mx=max(hopup[0].second,m_up[hopup[0].first]);
                for(size_t s=mn;s<mx;++s)
                    odd= (odd != (m_latupid[s]!=m_Nup));
                sup[m_up[hopup[0].first]]=m_Nup;
                sup[hopup[0].second]=hopup[0].first;
            }
            delete [] sup;
        }
    }
    if(hopdo.size()){
        mn=min(hopdo[0].second,m_do[hopdo[0].first])+1;
        mx=max(hopdo[0].second,m_do[hopdo[0].first]);
        for(size_t s=mn;s<mx;++s)
            odd= (odd != (m_latdoid[s]!=m_Ndo));
        if(hopdo.size()>1){
            size_t* sdo=new size_t[m_L*m_L];
            memcpy(sdo,m_latdoid,m_L*m_L*sizeof(size_t));
            sdo[m_do[hopdo[0].first]]=m_Ndo;
            sdo[hopdo[0].second]=hopdo[0].first;
            for(size_t h=1;h<hopdo.size();++h){
                mn=min(hopdo[h].second,m_do[hopdo[h].first])+1;
                mx=max(hopdo[h].second,m_do[hopdo[h].first]);
                for(size_t s=mn;s<mx;++s)
                    odd= (odd != (m_latdoid[s]!=m_Ndo));
                sdo[m_do[hopdo[h].first]]=m_Ndo;
                sdo[hopdo[h].second]=hopdo[h].first;
            }
            delete [] sdo;
        }
    }
    return (2*(!odd)-1);
}

int SpinState::cdagsign(size_t rr,bool up)
{
    int out=1;
    for(size_t r=0;r<rr;++r)
        if(m_latupid[r]!=m_Nup) out*=-1;
    if(!up){
        for(size_t r=0;r<rr;++r)
            if(m_latdoid[r]!=m_Ndo) out*=-1;
    }
    return out;
}

size_t SpinState::count_nn() const
{
    int taus[8]={1,0,-1,0,0,1,0,-1};
    int out(0);
    for(size_t x=0;x<m_L;++x){
        for(size_t y=0;y<m_L;++y){
            for(size_t t=0;t<4;++t){
                size_t xp,yp;
                xp=linalg::mod(x+taus[2*t],m_L);
                yp=linalg::mod(y+taus[2*t+1],m_L);
                out+=(m_latoc[y*m_L+x]!=m_latoc[yp*m_L+xp]);
            }
        }
    }
    return out/2;
}

std::ostream& operator<<(std::ostream& out,const SpinState& sp)
{
    out<<"********** SpinState: *************"<<std::endl;
    for(size_t y=0;y<sp.m_L;++y){
        for(size_t x=0;x<sp.m_L;++x){
            if(sp.m_latoc[y*sp.m_L+x]==EMPTY) out<<std::setw(2)<<" .";
            else if(sp.m_latoc[y*sp.m_L+x]==UP) out<<std::setw(2)<<"1";
            else if(sp.m_latoc[y*sp.m_L+x]==DOWN) out<<std::setw(2)<<"0";
            else out<<std::setw(2)<<"D";
        }
        out<<"   | ";
        for(size_t x=0;x<sp.m_L;++x){
            if(sp.m_latupid[y*sp.m_L+x]==sp.m_Nup)
                out<<"  .";
            else
                out<<" "<<std::setw(2)<<sp.m_latupid[y*sp.m_L+x];
        }
        out<<"   | ";
        for(size_t x=0;x<sp.m_L;++x){
            if(sp.m_latdoid[y*sp.m_L+x]==sp.m_Ndo)
                out<<"  .";
            else
                out<<" "<<std::setw(2)<<sp.m_latdoid[y*sp.m_L+x];
        }
        out<<std::endl;
    }
    out<<std::endl<<"Spins up:"<<std::endl;
    for(size_t up=0;up<sp.m_Nup;++up)
        out<<std::setw(2)<<sp.m_up[up]<<" ";
    out<<std::endl;
    out<<std::endl<<std::endl<<"Spins down:"<<std::endl;
    for(size_t down=0;down<sp.m_Ndo;++down)
        out<<std::setw(2)<<sp.m_do[down]<<" ";
    out<<std::endl<<std::endl;
    out<<"Sign: "<<sp.GetSign()<<std::endl;
    //out<<"upside-down NN: "<<sp.m_nn<<std::endl;
    out<<"***********************************"<<std::endl;
    return out;
}
