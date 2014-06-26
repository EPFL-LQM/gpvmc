#ifndef _JASTROWPOTENTIAL_H
#define _JASTROWPOTENTIAL_H

#include "defs.h"

class Lattice;

class JastrowPotential{
    public:
        JastrowPotential(const Lattice* lattice, const std::vector<double>& params);
        ~JastrowPotential();
        std::size_t GetNParams() const;
        double Pot(const uint_vec_t& statei,
                   const uint_vec_t& statej) const;
        void PotGrad(const uint_vec_t& statei,
                     const uint_vec_t& statej,
                     std::vector<double>& grad) const;
        void PotHess(const uint_vec_t& statei,
                     const uint_vec_t& statej,
                     std::vector<double>& hess) const;
    protected:
        virtual void init();
        virtual double space_potential(const uint_vec_t& Ri,
                                       const std::vector<double>& ri,
                                       const uint_vec_t& Rj,
                                       const std::vector<double>& rj) const=0;
        virtual double space_potential_grad(const uint_vec_t& Ri,
                                            const std::vector<double>& ri,
                                            const uint_vec_t& Rj,
                                            const std::vector<double>& rj,
                                            std::size_t param) const=0;
        virtual double space_potential_hess(const uint_vec_t& Ri,
                                            const std::vector<double>& ri,
                                            const uint_vec_t& Rj,
                                            const std::vector<double>& rj,
                                            std::size_t param_a,
                                            std::size_t param_b) const=0;
        virtual double internal_quantum_number_potential(const uint_vec_t& statei,
                                                         const uint_vec_t& statej) const=0;
        const Lattice* m_lattice;
        std::vector<double> m_params;
    private:
        double* m_rijpot;//!< Cache of the real space potential.
        std::vector<double*> m_rijpot_grad;//!< Cache of the real space potential gradient.
        std::vector<double*> m_rijpot_hess;//!< Cache of the real space potential hessian matrix.
};

#endif//_JASTROWPOTENTIAL_H
