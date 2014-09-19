#include "JastrowPotential.h"
#include "Lattice.h"
#include <cmath>

using namespace std;

JastrowPotential::JastrowPotential(const Lattice* lattice, const vector<double>& params)
    :m_lattice(lattice),
     m_params(params)
{
    m_rijpot=new double[lattice->GetNv()*lattice->GetNv()];
}

JastrowPotential::~JastrowPotential()
{
    delete [] m_rijpot;
}

void JastrowPotential::init()
{
    for(size_t vi=0;vi<m_lattice->GetNv();++vi){
        for(size_t vj=0;vj<m_lattice->GetNv();++vj){
            m_rijpot[vi*m_lattice->GetNv()+vj]=this->space_potential(
                    m_lattice->GetVertices()[vi]->uc,
                    m_lattice->GetVertices()[vi]->pos,
                    m_lattice->GetVertices()[vj]->uc,
                    m_lattice->GetVertices()[vj]->pos);
        }
    }
}

double JastrowPotential::Pot(const uint_vec_t& statei,const uint_vec_t& statej) const
{
    double out=m_rijpot[statei[0]*m_lattice->GetNv()+statej[0]]*
            this->internal_quantum_number_potential(statei,statej);
    return out;
}

void JastrowPotential::transvecmod(int & rx, int & ry) const
{
    rx=rx-(rx>=int(m_lattice->GetLx()/2))*int(m_lattice->GetLx())
         +(rx<-int(m_lattice->GetLx()/2))*int(m_lattice->GetLx());
    ry=ry-(ry>=int(m_lattice->GetLy()/2))*int(m_lattice->GetLy())
         +(ry<-int(m_lattice->GetLy()/2))*int(m_lattice->GetLy());
}
