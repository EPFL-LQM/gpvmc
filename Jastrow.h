#ifndef _JASTROW_H
#define _JASTROW_H

#include "defs.h"
#include "Amplitude.h"

class JastrowPotential;

class Jastrow: public Amplitude {
    public:
        Jastrow(const LatticeState* st,const JastrowPotential* pot);
        virtual void Init();
        virtual std::size_t GetNParams() const;
        virtual double Jas() const;
        virtual void JasGrad(std::vector<double>& grad) const;
        virtual void JasHess(std::vector<double>& hess) const;

        virtual void VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                                std::vector<double>& js) const;
        virtual void VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                                std::vector<double>& js,
                                std::vector<double>& js_grad,
                                std::vector<double>& js_hess) const;
        virtual void Update(const std::vector<hop_path_t>& rhop);
    protected:
        double m_jassum;
        std::vector<double> m_jassum_grad;
        std::vector<double> m_jassum_hess;
        const JastrowPotential* m_jaspot;
    private:
        void virtual_update(const std::vector<std::vector<hop_path_t> >& rhop,
                            std::vector<double>& js,
                            std::vector<double>& js_grad,
                            std::vector<double>& js_hess,
                            bool grad_hess) const;
};

#endif//_JASTROW_H
