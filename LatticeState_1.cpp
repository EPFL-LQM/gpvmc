#include "LatticeState_1.h"

using namespace std;

LatticeState_1::LatticeState_1(const Lattice* lattice,
                               const uint_vec_t& Npt,
                               const uint_vec_t& Nifs)
    : m_lattice(lattice), m_Nifs(Nifs)
{
    vector<size_t> lNfs(Npt.size(),m_lattice->vertices.size());
    for(size_t f=0;f<lNfs.size();++f) lNfs[f]*=Nifs[f];
    build_base(Npt,lNfs);
}

const Lattice* LatticeState_1::GetLattice() const
{
    return m_lattice;
}

void LatticeState_1::GetLatOc(size_t v,
                              vector<uint_vec_t>& st)
{
    st=vector<uint_vec_t>(m_Nfl);
    for(size_t fl=0;fl<m_Nfl;++fl){
        for(size_t i=0;i<m_Nifs[fl];++i){
            if(m_fock[fl][v*m_Nifs[fl]+i]!=m_Npt[fl]){
                st[fl].push_back(i);
            }
        }
    }
}

