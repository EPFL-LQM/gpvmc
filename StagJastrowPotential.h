#ifndef _STAGJASTROWPOTENTIAL_H
#define _STAGJASTROWPOTENTIAL_H

#include <iostream>
#include <algorithm>
#include "JastrowPotential.h"
#include "linalg.h"

using namespace std;

class StagJastrowPotential: public JastrowPotential {
    public:
        StagJastrowPotential(const Lattice* lattice, double neel)
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
            if(Ri[0]!=Rj[0] || Ri[1]!=Rj[1] || ri[0]!=rj[0] || ri[1]!=rj[1])
                out=m_params[0]*double(1.0-2.0*linalg::mod(int(Rj[0]+Rj[1])-int(Ri[0]+Ri[1]),2))/m_lattice->GetNv();
            return out;
        }

        virtual double internal_quantum_number_potential(const uint_vec_t& statei,
                                                         const uint_vec_t& statej) const
        {
            if(isup(statei)){
                if(isup(statej)){
                    return 1;
                } else {
                    return -1;
                }
            } else {
                if(isup(statej)){
                    return -1;
                } else {
                    return 1;
                }
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

#endif//_STAGJASTROWPOTENTIAL_H
