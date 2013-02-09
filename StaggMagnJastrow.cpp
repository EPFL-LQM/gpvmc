#include "StaggMagnJastrow.h"
#include "linalg.h"
#include "SpinState.h"
#include <iostream>

using namespace std;

StaggMagnJastrow::StaggMagnJastrow(SpinState* sp, double gamma)
    :Jastrow(sp), m_gamma(gamma)
{
    sp->SetJas(this);
    this->init();
}

void StaggMagnJastrow::init()
{
    m_sum=0;
    for(size_t xi=0;xi<m_sp->GetL();++xi){
        for(size_t yi=0;yi<m_sp->GetL();++yi){
            int sigi=2*int(m_sp->GetLatOc(xi,yi))-1;
            if(linalg::mod(xi+yi,2))
                m_sum+=sigi;
            else
                m_sum-=sigi;
        }
    }
    m_jas=std::exp(m_gamma*m_sum);
}

double StaggMagnJastrow::virtualhop(const hop_path_t& hopup,
                                    const hop_path_t& hopdo) const
{
    // Swap detection, other cases are not supported yet.
    double out(0);
    if(hopup.size()==1 && hopdo.size()==1 &&
       hopup[0].second==m_sp->GetDo()[hopdo[0].first] &&
       hopdo[0].second==m_sp->GetUp()[hopup[0].first])
    {
        size_t rup[2];
        size_t rdo[2];
        rup[1]=m_sp->GetUp()[hopup[0].first]/m_sp->GetL();
        rup[0]=m_sp->GetUp()[hopup[0].first]-rup[1]*m_sp->GetL();
        rdo[1]=m_sp->GetDo()[hopdo[0].first]/m_sp->GetL();
        rdo[0]=m_sp->GetDo()[hopdo[0].first]-rdo[1]*m_sp->GetL();
        out=4*(double(linalg::mod(rdo[0]+rdo[1],2))
                -double(linalg::mod(rup[0]+rup[1],2)));
        //return m_jas*std::exp(m_gamma*out/(m_sp->GetL()*m_sp->GetL()));
    } else {
#ifdef USE_EXC
        throw(logic_error("StaggMagnJastrow::virtualhop: "
                          "General hop not implemented yet "
                          "only up-down spin swap supported"));
#else
        cerr<<"StaggMagnJastrow::virtualhop: "
              "General hop not implemented yet "
              "only up-down spin swap supported"
            <<endl;
        abort();
#endif
    }
    return m_jas*std::exp(m_gamma*out);
    
}
