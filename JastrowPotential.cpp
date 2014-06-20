#include "JastrowPotential.h"
#include "Lattice.h"

using namespace std;

JastrowPotential::JastrowPotential(const Lattice* lattice)
    :m_lattice(lattice)
{
    m_rijpot=new double[lattice->GetNv()*lattice->GetNv()];
}

JastrowPotential::~JastrowPotential()
{
    delete [] m_rijpot;
}

void JastrowPotential::Init()
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
