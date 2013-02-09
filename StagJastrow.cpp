#include "StagJastrow.h"
#include "SpinState.h"
#include "linalg.h"

double StagJastrow::jpot(const size_t* r)
{
    double out(0);
    if(linalg::mod(r[0]+r[1],2))
        out=-1;
    else
        out=1;
    return m_gamma*out/std::pow(m_sp->GetL(),2);
}
