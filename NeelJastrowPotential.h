#ifndef _NEELJASTROWPOTENTIAL_H
#define _NEELJASTROWPOTENTIAL_H

#include <iostream>
#include <algorithm>
#include "JastrowPotential.h"

using namespace std;

class NeelJastrowPotential: public JastrowPotential {
    public:
        NeelJastrowPotential(const Lattice* lattice, double neel)
            :JastrowPotential(lattice,vector<double>(1,neel))
        {
            this->init();
        }
    protected:
        virtual double space_potential(const uint_vec_t& Ri,
                                       const std::vector<double>& ri,
                                       const uint_vec_t& Rj,
                                       const std::vector<double>& rj) const
        {
            double out(0);
            if(Ri[0]==Rj[0] && Ri[1]==Rj[1])
                out=double(m_params[0]*(1.0-2.0*(int(Ri[0]+Ri[1])%2)));
                //out=double(m_neel*(1.0-2.0*(int(Ri[0]+Ri[1])%2)))/sqrt(double(m_lattice->GetNv()));
            return out;
        }
        virtual double internal_quantum_number_potential(const uint_vec_t& statei,
                                                         const uint_vec_t& statej) const
        {
            if(equal(statei.begin(),statei.end(),statej.begin())){
                if(isup(statei))
                    return 1;
                else
                    return -1;
            } else {
                return 0;
            }
        }
    private:
        bool isup(const uint_vec_t& state) const
        {
            if(state[1]==0 and state[2]==0)
                return true;
            else
                return false;
        }
};

#endif//_NEELJASTROWPOTENTIAL_H
