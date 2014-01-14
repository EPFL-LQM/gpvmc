#include "SquareLattice.h"
#include "linalg.h"
#include <iomanip>
#include <sstream>
#include <initializer_list>

using namespace std;

SquareLattice::SquareLattice(size_t Lx, size_t Ly)
    : m_Lx(Lx), m_Ly(Ly)
{
    // By convention we store the lattice sites raw-wise
    // create vertices
    vector<Vertex*> non_const;
    for(size_t y=0;y<Ly;++y){
        for(size_t x=0;x<Lx;++x){
            non_const.push_back(new Vertex({ double(x), double(y) }));
            vertices.push_back(non_const.back());
        }
    }
    // create edges
    for(size_t x=0; x<Lx; ++x){
        for(size_t y=0; y<Ly; ++y){
            edges.push_back(new Edge(non_const[y*Lx+x],non_const[linalg::mod(y+1,Ly)*Lx+x]));
            edges.push_back(new Edge(non_const[y*Lx+x],non_const[y*Lx+linalg::mod(x+1,Lx)]));
        }
    }
}

string SquareLattice::str(std::vector<string> st) const
{
    if(!st.size())
        st=vector<string>(vertices.size());
    ostringstream out;
    for(size_t y=0;y<m_Ly;++y){
        for(size_t x=0;x<m_Lx;++x){
            out<<setfill('#')<<setw(3)<<st[y*m_Lx+x];
            if(x!=m_Lx-1) out<<'-';
        }
        if(y!=m_Ly-1){
            out<<endl;
            for(size_t x=0;x<m_Lx;++x){
                out<<" | ";
                if(x!=m_Lx-1) out<<" ";
            }
            out<<endl;
        }
    }
    return out.str();
}
