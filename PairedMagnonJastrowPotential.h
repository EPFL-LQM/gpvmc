#ifndef _PAIREDMAGNONJASTROWPOTENTIAL_H
#define _PAIREDMAGNONJASTROWPOTENTIAL_H

#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include "JastrowPotential.h"
#include "linalg.h"

class PairedMagnonJastrowPotential: public JastrowPotential {
    public:
        PairedMagnonJastrowPotential(const Lattice* lattice,
                                     const vector<double>& neigh_ratios);
        ~PairedMagnonJastrowPotential();
    protected:
        virtual void init();
        virtual double space_potential(const uint_vec_t& Ri,
                                       const std::vector<double>& ri,
                                       const uint_vec_t& Rj,
                                       const std::vector<double>& rj) const;
        virtual double internal_quantum_number_potential(const uint_vec_t& statei,
                                                         const uint_vec_t& statej) const;
    private:
        bool isup(const uint_vec_t& state) const;

        map<size_t,size_t> m_neigh_index;
        double* m_cache;
        double* m_neigh_r;
};

#endif//_LOGJASTROWPOTENTIAL_H
