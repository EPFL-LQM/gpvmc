#ifndef _STAGGMAGNJASTROW_H
#define _STAGGMAGNJASTROW_H

#include "Jastrow.h"

class StaggMagnJastrow: public Jastrow {
    public:
        StaggMagnJastrow(SpinState* sp,double gamma);
        double Jas() const {return m_jas;}
        virtual double virtualhop(const hop_path_t& hopup,
                                  const hop_path_t& hopdo) const;
        virtual void init();
        friend std::ostream& operator<<(std::ostream& os,const Jastrow& j);
    protected:
        virtual double jpot(const size_t* r){return 0;}
        double m_sum;
        double m_gamma;
};

#endif//_STAGGMAGNJASTROW_H
