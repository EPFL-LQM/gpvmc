#ifndef _LOGJASTROWPOTENTIAL_H
#define _LOGJASTROWPOTENTIAL_H

#include <iostream>
#include <algorithm>
#include "JastrowPotential.h"
#include "linalg.h"

using namespace std;

class LogJastrowPotential: public JastrowPotential {
    public:
        LogJastrowPotential(const Lattice* lattice, double neel)
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
            if(Ri[0]==Rj[0] && Ri[1]==Rj[1])
                return 0.0;
            else{
                int rx,ry;
                rx=int(Rj[0])-int(Ri[0]);
                ry=int(Rj[1])-int(Ri[1]);
                int ph=1;
                if(linalg::mod(rx+ry,2))
                    ph=-1;
                return -ph*m_params[0]*log(sqrt(0.5*(pow(sin(rx*M_PI/m_lattice->GetLx()),2)+pow(sin(ry*M_PI/m_lattice->GetLy()),2))))/m_lattice->GetNv();
            }
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

#endif//_LOGJASTROWPOTENTIAL_H
