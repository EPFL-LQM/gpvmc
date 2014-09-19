#ifndef _JASTROWPOTENTIAL_H
#define _JASTROWPOTENTIAL_H

#include "defs.h"

class Lattice;

class JastrowPotential{
    public:
        JastrowPotential(const Lattice* lattice, const std::vector<double>& params);
        ~JastrowPotential();
        double Pot(const uint_vec_t& statei,
                   const uint_vec_t& statej) const;
    protected:
        virtual void init();
        virtual double space_potential(const uint_vec_t& Ri,
                                       const std::vector<double>& ri,
                                       const uint_vec_t& Rj,
                                       const std::vector<double>& rj) const=0;
        virtual double internal_quantum_number_potential(const uint_vec_t& statei,
                                                         const uint_vec_t& statej) const=0;
        void transvecmod(int& rx, int& ry) const;
        const Lattice* m_lattice;
        std::vector<double> m_params;
    private:
        double* m_rijpot;//!< Cache of the real space potential.
};

#endif//_JASTROWPOTENTIAL_H
