#ifndef _PAIREDMAGNONJASTROWPOTENTIAL_H
#define _PAIREDMAGNONJASTROWPOTENTIAL_H

#include <iostream>
#include <algorithm>
#include "JastrowPotential.h"
#include "linalg.h"

using namespace std;

class PairedMagnonJastrowPotential: public JastrowPotential {
    public:
        PairedMagnonJastrowPotential(const Lattice* lattice, double neel, double nnjas, double nnnjas)
            :JastrowPotential(lattice,vector<double>({neel,nnjas,nnnjas}))
        {
            this->init();
        }
    protected:
        virtual double space_potential(const uint_vec_t& Ri,
                                       const std::vector<double>& ri,
                                       const uint_vec_t& Rj,
                                       const std::vector<double>& rj) const
        {
            double rx,ry;
            rx=double(Rj[0])-double(Ri[0]);
            ry=double(Rj[1])-double(Ri[1]);
            if(Ri[0]==Rj[0] && Ri[1]==Rj[1]){
                return 0.0;
            } else if(NN(Ri,Rj)){
                return -m_params[1]*func(rx,ry);
            } else if(NNN(Ri,Rj)){
                return -m_params[2]*func(rx,ry);
            } else{
                return -m_params[0]*func(rx,ry);
            }
        }
        virtual double space_potential_grad(const uint_vec_t& Ri,
                                            const std::vector<double>& ri,
                                            const uint_vec_t& Rj,
                                            const std::vector<double>& rj,
                                            std::size_t param) const
        {
            double rx,ry;
            rx=double(Rj[0])-double(Ri[0]);
            ry=double(Rj[1])-double(Ri[1]);
            if(Ri[0]==Rj[0] && Ri[1]==Rj[1]){
                return 0.0;
            } else if(NN(Ri,Rj)){
                if(param==1) return -func(rx,ry);
                else return 0.0;
            } else if(NNN(Ri,Rj)){
                if(param==2) return -func(rx,ry);
                else return 0.0;
            } else{
                return -func(rx,ry);
            }
        }
        virtual double space_potential_hess(const uint_vec_t& Ri,
                                            const std::vector<double>& ri,
                                            const uint_vec_t& Rj,
                                            const std::vector<double>& rj,
                                            std::size_t param_a,
                                            std::size_t param_b) const
        {
            return 0;
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
        double func(const double& rx, const double& ry) const {
            // first determine the k=(0,0) contribution.
            // It is chosen such that func(m_lattice.Lx,m_lattice.Ly)=0
            // and corresponds to a global shift of the potential.
            double shift=0;
            for(size_t kx=0;kx<m_lattice->GetLx();++kx){
                for(size_t ky=0;ky<m_lattice->GetLy();++ky){
                    if(kx!=0 && ky!=0){
                        double dkx=double(kx)*2*M_PI/double(m_lattice->GetLx());
                        double dky=double(ky)*2*M_PI/double(m_lattice->GetLy());
                        double gk=0.5*(cos(dkx)+cos(dky));
                        shift+=(sqrt((1+gk)/(1-gk))-1)*cos(dkx*m_lattice->GetLx()+dky*m_lattice->GetLy());//don't include complex part since gk is even
                    }
                }
            }
            // calculate for rx and ry:
            double out=0;
            for(size_t kx=0;kx<m_lattice->GetLx();++kx){
                for(size_t ky=0;ky<m_lattice->GetLy();++ky){
                    if(kx!=0 && ky!=0){
                        double dkx=double(kx)*2*M_PI/double(m_lattice->GetLx());
                        double dky=double(ky)*2*M_PI/double(m_lattice->GetLy());
                        double gk=0.5*(cos(dkx)+cos(dky));
                        out+=(sqrt((1+gk)/(1-gk))-1)*cos(dkx*rx+dky*ry);//don't include complex part since gk is even
                    }
                }
            }
            return (out-shift)/m_lattice->GetLx()/m_lattice->GetLy();
        }
        bool NN(const uint_vec_t& Ri,const uint_vec_t& Rj) const
        {
            return (abs(int(Ri[0])-int(Rj[0]))==1 && Ri[1]==Rj[1]) ||
                   (abs(int(Ri[1])-int(Rj[1]))==1 && Ri[0]==Rj[0]);
        }
        bool NNN(const uint_vec_t& Ri,const uint_vec_t& Rj) const
        {
            return abs(int(Ri[0])-int(Rj[0]))==1 &&
                   abs(int(Ri[1])-int(Rj[1]))==1;
        }
};

#endif//_LOGJASTROWPOTENTIAL_H
