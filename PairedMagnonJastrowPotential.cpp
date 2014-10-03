#include "PairedMagnonJastrowPotential.h"
#include "Lattice.h"
#include <fftw3.h>
#include <fstream>

using namespace std;

PairedMagnonJastrowPotential::PairedMagnonJastrowPotential(
        const Lattice* lattice,
        const vector<double>& neigh_ratio)
    :JastrowPotential(lattice,neigh_ratio),
     m_cache(new double[(lattice->GetLx()/2+1)*(lattice->GetLy()/2+1)])
{
    this->init();
    std::ofstream fdb;
    fdb.open("pa_pot_db.txt");
    for(int x=-int(m_lattice->GetLx()/2);x<int(m_lattice->GetLx()/2);++x){
        for(int y=-int(m_lattice->GetLy()/2);y<int(m_lattice->GetLy()/2);++y){
            fdb<<m_cache[abs(x)*(m_lattice->GetLx()/2+1)+abs(y)]<<"\t";
        }
        fdb<<std::endl;
    }
    fdb.close();
}

PairedMagnonJastrowPotential::~PairedMagnonJastrowPotential()
{
    delete [] m_cache;
    delete [] m_neigh_r;
}

void PairedMagnonJastrowPotential::init()
{
    size_t Lx(m_lattice->GetLx()),Ly(m_lattice->GetLy());
    // Initialize the neighbour index set.
    set<size_t> neigh_dist;
    for(size_t rx=0;rx<=Lx/2;++rx){
        for(size_t ry=0;ry<=Ly/2;++ry){
            neigh_dist.insert(rx*rx+ry*ry);
        }
    }
    size_t idx=0;
    for(set<size_t>::iterator iter=neigh_dist.begin();iter!=neigh_dist.end();++iter){
        if(idx<m_params.size()){
            m_neigh_index[*iter]=idx;
            idx++;
        } else {
            m_neigh_index[*iter]=m_params.size()-1;
        }
    }

    double * Dk = new double[(Lx/2+1)*(Ly/2+1)];
    fftw_plan plan = fftw_plan_r2r_2d(Lx/2+1,Ly/2+1,
                                      Dk,m_cache,
                                      FFTW_REDFT00,FFTW_REDFT00,
                                      FFTW_ESTIMATE);
    Dk[0]=0;
    for(size_t kx=0;kx<=Lx/2;++kx){
        for(size_t ky=0;ky<=Ly/2;++ky){
            if(kx!=0 || ky!=0){
                double dkx=double(kx)*2*M_PI/double(Lx);
                double dky=double(ky)*2*M_PI/double(Ly);
                double gk=0.5*(cos(dkx)+cos(dky));
                Dk[kx*(Ly/2+1)+ky]=sqrt((1+gk)/(1-gk))-1;
            }
        }
    }
    fftw_execute(plan);
    for(size_t x=0;x<=Lx/2;++x){
        for(size_t y=0;y<=Ly/2;++y){
            m_cache[x*(Ly/2+1)+y]-=m_cache[(Lx/2)*(Ly/2+1)+Ly/2];
            m_cache[x*(Ly/2+1)+y]/=Lx*Ly;
        }
    }
    fftw_destroy_plan(plan);
    delete [] Dk;
    JastrowPotential::init();
}

double PairedMagnonJastrowPotential::space_potential(const uint_vec_t& Ri,
                                                     const std::vector<double>& ri,
                                                     const uint_vec_t& Rj,
                                                     const std::vector<double>& rj
                                                     ) const
{
    int rx,ry;
    rx=int(Rj[0])-int(Ri[0]);
    ry=int(Rj[1])-int(Ri[1]);
    transvecmod(rx,ry);
    return -m_params[m_neigh_index.at(abs(rx)*(m_lattice->GetLy()/2+1)+abs(ry))]*
        m_cache[abs(rx)*(m_lattice->GetLy()/2+1)+abs(ry)];
}
double PairedMagnonJastrowPotential::internal_quantum_number_potential(
                                                const uint_vec_t& statei,
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
bool PairedMagnonJastrowPotential::isup(const uint_vec_t& state) const
{
    if(state[1]==0 and state[2]==0)
        return true;
    else
        return false;
}
