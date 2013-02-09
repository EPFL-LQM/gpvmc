#include "Jastrow.h"
#include "linalg.h"
#include "SpinState.h"

using namespace std;

Jastrow::Jastrow(SpinState* sp)
    :m_sp(sp),
     m_pot(new double[m_sp->GetL()*m_sp->GetL()]),
     m_jas(1)
{
    sp->SetJas(this);
}

Jastrow::~Jastrow()
{
    delete [] m_pot;
}

void Jastrow::init()
{
    if(m_pot) delete [] m_pot;
    m_pot=new double[m_sp->GetL()*m_sp->GetL()];
    size_t r[2];
    for(r[0]=0;r[0]<m_sp->GetL();++r[0]){
        for(r[1]=0;r[1]<m_sp->GetL();++r[1]){
            m_pot[r[0]*m_sp->GetL()+r[1]]=this->jpot(r);
        }
    }
    double out(0);
    for(size_t xi=0;xi<m_sp->GetL();++xi){
        for(size_t yi=0;yi<m_sp->GetL();++yi){
            for(size_t xj=0;xj<m_sp->GetL();++xj){
                for(size_t yj=0;yj<m_sp->GetL();++yj){
                    size_t dx=linalg::mod(int(xi)-int(xj),m_sp->GetL());
                    size_t dy=linalg::mod(int(yi)-int(yj),m_sp->GetL());
                    int sigi=2*int(m_sp->GetLatOc(xi,yi))-1;
                    int sigj=2*int(m_sp->GetLatOc(xj,yj))-1;
                    out+=m_pot[dx*m_sp->GetL()+dy]*sigi*sigj;
                }
            }
        }
    }
    m_jas=exp(0.5*out);
}

double Jastrow::virtualhop(const hop_path_t& hopup,
                           const hop_path_t& hopdo) const
{
    double out(0);
    vector<size_t> iup, ido;
    for(size_t h=0;h<hopup.size();++h){
        iup.push_back(hopup[h].second);
        iup.push_back(m_sp->GetUp()[hopup[h].first]);
    }
    for(size_t h=0;h<hopdo.size();++h){
        ido.push_back(hopdo[h].second);
        ido.push_back(m_sp->GetDo()[hopdo[h].first]);
    }
    // swapping detection for simplest update
    if(hopup.size()==1 && hopdo.size()==1 &&
            iup[0]==ido[1] && iup[1]==ido[0]){
        size_t drup[2], drdo[2], rup[2], rdo[2];
        rup[0]=iup[1]%m_sp->GetL();
        rup[1]=iup[1]/m_sp->GetL();
        rdo[0]=iup[0]%m_sp->GetL();
        rdo[1]=iup[0]/m_sp->GetL();
        for(size_t xi=0;xi<m_sp->GetL();++xi){
            for(size_t yi=0;yi<m_sp->GetL();++yi){
                if((xi!=rup[0] || yi!=rup[1]) &&
                   (xi!=rdo[0] || yi!=rdo[1])){
                    drup[0]=linalg::mod(int(xi)-int(rup[0]),
                                        m_sp->GetL());
                    drup[1]=linalg::mod(int(yi)-int(rup[1]),
                                        m_sp->GetL());
                    drdo[0]=linalg::mod(int(xi)-int(rdo[0]),
                                        m_sp->GetL());
                    drdo[1]=linalg::mod(int(yi)-int(rdo[1]),
                                        m_sp->GetL());
                    out+=(m_pot[drdo[0]*m_sp->GetL()+drdo[1]]-
                          m_pot[drup[0]*m_sp->GetL()+drup[1]])*
                          (2*int(m_sp->GetLatOc(xi,yi))-1);
                }
            }
        }
    } else {
#ifdef USE_EXC
        throw logic_error("Jastrow::virtualhop: General "
                          "hop sequence not implemented, "
                          "only swapping works...");
#else
        cerr<<"Jastrow::virtualhop: General "
              "hop sequence not implemented, "
              "only swapping works..."<<endl;
        abort();
#endif
    }
    return m_jas*exp(2*out);
}

void Jastrow::hop(const hop_path_t& hopup,
                  const hop_path_t& hopdo)
{
    m_jas=virtualhop(hopup,hopdo);
}

std::ostream& operator<<(std::ostream& os,const Jastrow& j)
{

    os<<std::setw(8)<<" "<<" ";
    for(size_t y=0;y<j.m_sp->GetL();++y)
        os<<std::setw(8)<<y<<" ";
    os<<std::endl;
    for(size_t x=0;x<j.m_sp->GetL();++x){
        os<<std::setw(8)<<x<<" ";
        for(size_t y=0;y<j.m_sp->GetL();++y)
            os<<std::setw(8)<<j.m_pot[x*j.m_sp->GetL()+y]<<" ";
        os<<std::endl;
    }
    return os;
}
