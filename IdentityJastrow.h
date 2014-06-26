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
        virtual std::size_t GetNParams() const
        {
            return 0;
        }
        virtual double Jas() const
        {
            return 1;
        }
        virtual void JasGrad(std::vector<double>& grad) const
        {
            grad.resize(0);
        }
        virtual void JasHess(std::vector<double>& hess) const
        {
            hess.resize(0);
        }
        virtual void VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                                std::vector<double>& js) const
        {
            js.assign(rhop.size(),1.0);
        }
        virtual void VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                                std::vector<double>& js,
                                std::vector<double>& js_grad,
                                std::vector<double>& js_hess) const
        {
            js.assign(rhop.size(),1.0);
            js_grad.resize(0);
            js_hess.resize(0);
        }
        virtual void Update(const std::vector<hop_path_t>& rhop)
        {}
};

#endif//_IDENTITYJASTROW_H
