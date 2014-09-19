#ifndef _PAIREDMAGNONJASTROWPOTENTIAL_H
#define _PAIREDMAGNONJASTROWPOTENTIAL_H

#include <iostream>
#include <algorithm>
#include "JastrowPotential.h"
#include "linalg.h"

class PairedMagnonJastrowPotential: public JastrowPotential {
    public:
        PairedMagnonJastrowPotential(const Lattice* lattice,
                                     double neel, double nnjas,
                                     double nnnjas);
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
        bool NN(const double & rx,const double& ry) const;
        bool NNN(const double & rx,const double& ry) const;

        double* m_cache;
};

#endif//_LOGJASTROWPOTENTIAL_H
