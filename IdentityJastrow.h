#ifndef _IDENTITYJASTROW_H
#define _IDENTITYJASTROW_H

#include "Jastrow.h"

class IdentityJastrow: public Jastrow {
    public:
        IdentityJastrow()
            :Jastrow(0,0) {}
        virtual void Init()
        {
            m_jassum=0;
        }
        virtual void VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                                std::vector<double>& js) const
        {
            js.assign(rhop.size(),1.0);
        }
        virtual void Update(const std::vector<hop_path_t>& rhop)
        {}
};

#endif//_IDENTITYJASTROW_H
