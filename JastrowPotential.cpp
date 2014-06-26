#include "JastrowPotential.h"
#include "Lattice.h"
#include <cmath>

using namespace std;

JastrowPotential::JastrowPotential(const Lattice* lattice, const vector<double>& params)
    :m_lattice(lattice),
     m_params(params),
     m_rijpot_grad(params.size()),
     m_rijpot_hess(pow(params.size(),2))
{
    m_rijpot=new double[lattice->GetNv()*lattice->GetNv()];
    for(size_t pa=0;pa<m_params.size();++pa){
        m_rijpot_grad[pa]=new double[lattice->GetNv()*lattice->GetNv()];
        for(size_t pb=0;pb<m_params.size();++pb)
            m_rijpot_hess[pa*m_params.size()+pb]=new double[lattice->GetNv()*lattice->GetNv()];
    }
}

JastrowPotential::~JastrowPotential()
{
    delete [] m_rijpot;
    for(size_t pa=0;pa<m_params.size();++pa){
        delete [] m_rijpot_grad[pa];
        for(size_t pb=0;pb<m_params.size();++pb)
            delete [] m_rijpot_hess[pa*m_params.size()+pb];
    }
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
            for(size_t pa=0;pa<m_params.size();++pa){
                m_rijpot_grad[pa][vi*m_lattice->GetNv()+vj]=this->space_potential_grad(
                    m_lattice->GetVertices()[vi]->uc,
                    m_lattice->GetVertices()[vi]->pos,
                    m_lattice->GetVertices()[vj]->uc,
                    m_lattice->GetVertices()[vj]->pos,
                    pa);
                for(size_t pb=0;pb<m_params.size();++pb){
                    m_rijpot_hess[pa*m_params.size()+pb][vi*m_lattice->GetNv()+vj]=
                        this->space_potential_hess(
                        m_lattice->GetVertices()[vi]->uc,
                        m_lattice->GetVertices()[vi]->pos,
                        m_lattice->GetVertices()[vj]->uc,
                        m_lattice->GetVertices()[vj]->pos,
                        pa,pb);
                }
            }
        }
    }
}

size_t JastrowPotential::GetNParams() const
{
    return m_params.size();
}

double JastrowPotential::Pot(const uint_vec_t& statei,const uint_vec_t& statej) const
{
    double out=m_rijpot[statei[0]*m_lattice->GetNv()+statej[0]]*
            this->internal_quantum_number_potential(statei,statej);
    return out;
}

void JastrowPotential::PotGrad(const uint_vec_t& statei,
                               const uint_vec_t& statej,
                               vector<double>& grad) const
{
    grad.resize(m_params.size());
    for(size_t p=0;p<m_params.size();++p)
        grad[p]=m_rijpot_grad[p][statei[0]*m_lattice->GetNv()+statej[0]]*
            this->internal_quantum_number_potential(statei,statej);
}

void JastrowPotential::PotHess(const uint_vec_t& statei,
                               const uint_vec_t& statej,
                               vector<double>& hess) const
{
    hess.resize(pow(m_params.size(),2));
    for(size_t pa=0;pa<m_params.size();++pa)
        for(size_t pb=0;pb<m_params.size();++pb)
            hess[pa*m_params.size()+pb]=
                m_rijpot_hess[pa*m_params.size()+pb][statei[0]*m_lattice->GetNv()+statej[0]]*
                this->internal_quantum_number_potential(statei,statej);
}
