#include "LatticeState.h"
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>
#include "RanGen.h"
#include "Lattice.h"

using namespace std;

LatticeState::LatticeState(FileManager* fm,
                           const Lattice* lattice,
                           const uint_vec_t& Npt,
                           const uint_vec_t& Nifs)
    : State(fm), m_lattice(lattice), m_Nifs(Nifs)
{
    vector<size_t> lNfs(Npt.size(),m_lattice->GetNv());
    for(size_t f=0;f<lNfs.size();++f) lNfs[f]*=Nifs[f];
    build_base(Npt,lNfs);
}

size_t LatticeState::GetNsites() const
{
    return m_lattice->GetNv();
}

const uint_vec_t& LatticeState::GetNifs() const
{
    return m_Nifs;
}

const Lattice* LatticeState::GetLattice() const
{
    return m_lattice;
}

void LatticeState::RanInit(const vector<vector<size_t> >& pop)
{
    vector<uint_vec_t> fst(m_Nfl);
    for(size_t f=0; f<m_Nfl;++f){
        fst[f]=uint_vec_t(m_Nfs[f],m_Npt[f]);
    }
    for(size_t f=0; f<m_Nfl;++f){
        size_t np=0;
        for(size_t ni=0;ni<m_Nifs[f];++ni){
            size_t NTpt=0;
            for_each(m_Npt.begin(),m_Npt.end(),[&](size_t n){NTpt+=n;});
            if(f==m_Nfl-1 && ni==m_Nifs[f]-1 && NTpt==m_lattice->GetNv()){
                for(size_t v=0;v<m_lattice->GetNv();++v){
                    bool occupied=false;
                    for(size_t fp=0;fp<m_Nfl;++fp){
                        for(size_t si=0;si<m_Nifs[fp];++si){
                            if(fst[fp][v*m_Nifs[fp]+si]<m_Npt[fp])
                                occupied=true;
                        }
                    }
                    if(!occupied){
                        fst[f][v*m_Nifs[f]+ni]=np;
                        ++np;
                    }
                }
            } else {
                while(np!=pop[f][ni]){
                    size_t v=size_t(RanGen::uniform()*m_lattice->GetNv());
                    // check whether any particle occupies site v.
                    bool occupied=false;
                    for(size_t fp=0;fp<m_Nfl;++fp){
                        for(size_t si=0;si<m_Nifs[fp];++si){
                            if(fst[fp][v*m_Nifs[fp]+si]<m_Npt[fp])
                                occupied=true;
                        }
                    }
                    if(!occupied){
                        fst[f][v*m_Nifs[f]+ni]=np;
                        ++np;
                    }
                }
            }
        }
    }
    InitFock(fst);
}

void LatticeState::FockInit(const vector<vector<int> >& fock)
{
    bool compatible=true;
    if(fock.size()!=m_Nfl)
        compatible=false;
    else{
        for(size_t fl=0;fl<m_Nfl;++fl){
            if(m_Nfs[fl]!=fock[fl].size())
                compatible=false;
        }
    }
    if(!compatible){
#ifdef EXCEPT
        throw(std::runtime_error("LatticeState::FockInit: Fock state has no matching"
                                 " nuumber of flavours or different sizes"));
#else
        cerr<<"LatticeState::FockInit: Fock state has no matching"
              " nuumber of flavours or different sizes"<<endl;
        abort();
#endif
    }
    vector<uint_vec_t> fst(fock.size());
    for(size_t fl=0;fl<m_Nfl;++fl){
        fst[fl]=vector<size_t>(m_Nfs[fl],m_Npt[fl]);
        for(size_t f=0;f<m_Nfs[fl];++f){
            if(fock[fl][f]) fst[fl][f]=0;
        }
    }
    InitFock(fst);
}

void LatticeState::GetLatOc(size_t v,
                              vector<uint_vec_t>& st) const
{
    st.clear();
    st.resize(m_Nfl);
    for(size_t fl=0;fl<m_Nfl;++fl){
        for(size_t i=0;i<m_Nifs[fl];++i){
            if(m_fock[fl][v*m_Nifs[fl]+i]!=m_Npt[fl]){
                st[fl].push_back(i);
            }
        }
    }
}

void LatticeState::Fock2QN(size_t fidx, size_t flav, uint_vec_t& st) const
{
    // a fock index is constructed as follow:
    // fidx=v*m_Nifs[flav]+is
    // where v is the vertex index and is is the vertex
    // internal degree of freedom (for instance spin)
    // This is restricted to 2D and Bravais lattice so far...
    st.resize(3);
    st[0]=fidx/m_Nifs[flav];
    st[1]=fidx%m_Nifs[flav];
    st[2]=flav;
}

ostream& operator<<(ostream& out,const LatticeState& lst){
    vector<string> st(lst.m_lattice->GetNv());
    for(size_t v=0;v<lst.m_lattice->GetNv();++v){
        ostringstream s;
        size_t idx(0);
        for(size_t fl=0;fl<lst.m_Nfl;++fl){
            for(size_t i=0;i<lst.m_Nifs[fl];++i){
                if(lst.m_fock[fl][v*lst.m_Nifs[fl]+i]!=lst.m_Npt[fl])
                    idx+=pow(2,fl*(*max_element(lst.m_Nifs.begin(),lst.m_Nifs.end()))+i);
            }
        }
        s<<idx;
        st[v]=s.str();
    }
    out<<(State&)lst<<endl<<endl<<"Lattice state:"<<endl<<endl<<lst.m_lattice->str(st)<<endl;
    return out;
}
