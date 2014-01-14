#include "LatticeState_1.h"
#include <sstream>
#include <string>
#include <algorithm>

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

size_t LatticeState_1::GetNsites() const
{
    return m_lattice->vertices.size();
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

ostream& operator<<(ostream& out,const LatticeState_1& lst){
    vector<string> st(lst.m_lattice->vertices.size());
    for(size_t v=0;v<lst.m_lattice->vertices.size();++v){
        ostringstream s;
        size_t idx(0);
        for(size_t fl=0;fl<lst.m_Nfl;++fl){
            for(size_t i=0;i<lst.m_Nifs[fl];++i){
                if(lst.m_fock[fl][v*lst.m_Nifs[fl]]!=lst.m_Npt[fl])
                    idx+=pow(2,fl*(*max_element(lst.m_Nifs.begin(),lst.m_Nifs.end()))+i);
            }
        }
        s<<idx;
        st[v]=s.str();
    }
    out<<(State_1&)lst<<endl<<endl<<"Lattice state:"<<endl<<endl<<lst.m_lattice->str(st)<<endl;
    return out;
}
