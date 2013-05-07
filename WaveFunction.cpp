#include "WaveFunction.h"
#include "Amplitude.h"
#include "FileManager.h"
#include <hdf5_hl.h>
#include <hdf5.h>
#include <cstring>
#ifdef USEMPI
#include <mpi.h>
#endif

using namespace std;

WaveFunction::WaveFunction(size_t Lx, size_t Ly,
                           size_t Nfsup, size_t Nfsdo)
    : m_fsup(new size_t[Nfsup]),
      m_fsdo(new size_t[Nfsdo]),
      m_fockup(new size_t[Lx*Ly]),
      m_fockdo(new size_t[Lx*Ly]),
      m_cacheup(new std::complex<double>[Lx*Lx*Ly*Ly]),
      m_cachedo(new std::complex<double>[Lx*Lx*Ly*Ly]),
      m_Nfsup(Nfsup), m_Nfsdo(Nfsdo),
      m_Lx(Lx), m_Ly(Ly), m_sign(1)
{
    // Lx and Ly specify real space lattice
    // size which generate a complete
    // basis. The Fock space size is thus 2*Lx*Ly.
    for(size_t s=0;s<Lx*Ly;++s){
        // m_Nfsup means empty for spins up.
        m_fockup[s]=m_Nfsup;
        // m_Nfsdo means empty for spins down.
        m_fockdo[s]=m_Nfsdo;
    }
}

WaveFunction::~WaveFunction()
{
    delete [] m_fsup;
    delete [] m_fsdo;
    delete [] m_fockup;
    delete [] m_fockdo;
    delete [] m_cacheup;
    delete [] m_cachedo;
}

void WaveFunction::init_matrices()
{
    for(size_t f=0;f<m_Lx*m_Ly;++f){
        for(size_t R=0;R<m_Lx*m_Ly;++R){
            m_cacheup[f*m_Lx*m_Ly+R]=
                this->matrix_element(f,R,true);
            m_cachedo[f*m_Lx*m_Ly+R]=
                this->matrix_element(f,R,false);
        }
    }
}

int WaveFunction::hop_sign(const hop_path_t& hopup,
                           const hop_path_t& hopdo) const
{
    // calculate sign.
    bool odd=false;
    size_t mn,mx;
    if(hopup.size()){
        mn=min(hopup[0].second,m_fsup[hopup[0].first])+1;
        mx=max(hopup[0].second,m_fsup[hopup[0].first]);
        for(size_t s=mn;s<mx;++s)
            odd= (odd != (m_fockup[s]!=m_Nfsup));
        if(hopup.size()>1){
            size_t* sup=new size_t[m_Lx*m_Ly];
            std::memcpy(sup,m_fockup,m_Lx*m_Ly*sizeof(size_t));
            sup[m_fsup[hopup[0].first]]=m_Nfsup;
            sup[hopup[0].second]=hopup[0].first;
            for(size_t h=0;h<hopup.size();++h){
                mn=min(hopup[h].second,m_fsup[hopup[h].first])+1;
                mx=max(hopup[h].second,m_fsup[hopup[h].first]);
                for(size_t s=mn;s<mx;++s)
                    odd = (odd != (sup[s] != m_Nfsup));
                sup[m_fsup[hopup[h].first]]=m_Nfsup;
                sup[hopup[h].second]=hopup[h].first;
            }
            delete [] sup;
        }
    }
    if(hopdo.size()){
        mn=min(hopdo[0].second,m_fsdo[hopdo[0].first])+1;
        mx=max(hopdo[0].second,m_fsdo[hopdo[0].first]);
        for(size_t s=mn;s<mx;++s)
            odd= (odd != (m_fockdo[s]!=m_Nfsdo));
        if(hopdo.size()>1){
            size_t* sdo=new size_t[m_Lx*m_Ly];
            std::memcpy(sdo,m_fockdo,m_Lx*m_Ly*sizeof(size_t));
            sdo[m_fsdo[hopdo[0].first]]=m_Nfsdo;
            sdo[hopdo[0].second]=hopdo[0].first;
            for(size_t h=0;h<hopdo.size();++h){
                mn=min(hopdo[h].second,m_fsdo[hopdo[h].first])+1;
                mx=max(hopdo[h].second,m_fsdo[hopdo[h].first]);
                for(size_t s=mn;s<mx;++s)
                    odd = (odd != (sdo[s] != m_Nfsdo));
                sdo[m_fsdo[hopdo[h].first]]=m_Nfsdo;
                sdo[hopdo[h].second]=hopdo[h].first;
            }
            delete [] sdo;
        }
    }
    return (2*(!odd)-1);
}

void WaveFunction::add_state(const size_t* fup,
                             const size_t* fdo)
{
    if(!m_excup.size()){
        size_t stiup=0, stido=0;
        for(size_t s=0;s<m_Lx*m_Ly;++s){
            if(fup[s]!=m_Nfsup){
                m_fockup[s]=stiup;
                m_fsup[stiup]=s;
                ++stiup;
            }
            if(fdo[s]!=m_Nfsdo){
                m_fockdo[s]=stido;
                m_fsdo[stido]=s;
                ++stido;
            }
        }
        //arbitrarily assign +1 sign for the first state.
        m_sign=1;
        m_c_exc=0;
    }
    hop_path_t em;// empty hop
    // span hop path matrix with empty hops
    for(size_t i=0;i<m_excup.size();++i)
        m_excup[i].push_back(em);
    m_excup.push_back(
            vector<hop_path_t>(m_excup.size()+1,em));
    vector<size_t> state(m_Lx*m_Ly,0);
    for(size_t i=0;i<m_excup.size()-1;++i){
        // define line m_excup[end,i]
        get_hop_state(m_excup[0][i],m_fockup,m_Nfsup,state);
        m_excup.back()[i]=get_hop_path(
                fup,&state[0],m_Nfsup);
        // define line m_excup[i,end]
        m_excup[i].back()=get_hop_path(
                &state[0],fup,m_Nfsup);
    }
    // span matrix with empty hops
    for(size_t i=0;i<m_excdo.size();++i)
        m_excdo[i].push_back(em);
    vector<hop_path_t> temp(m_excdo.size()+1,em);
    m_excdo.push_back(temp);
    for(size_t i=0;i<m_excdo.size()-1;++i){
        // define line m_excdo[end,i]
        get_hop_state(m_excdo[0][i],m_fockdo,m_Nfsdo,state);
        m_excdo.back()[i]=get_hop_path(
                fdo,&state[0],m_Nfsdo);
        // define line m_excdo[i,end]
        m_excdo[i].back()=get_hop_path(
                &state[0],fdo,m_Nfsdo);
    }
    m_fock_states_up.push_back(vector<int>(m_Lx*m_Ly));
    m_fock_states_do.push_back(vector<int>(m_Lx*m_Ly));
    for(size_t f=0;f<m_Lx*m_Ly;++f){
        m_fock_states_up.back()[f]=(fup[f]!=m_Nfsup);
        m_fock_states_do.back()[f]=(fdo[f]!=m_Nfsdo);
    }
}

void WaveFunction::get_hop_state(hop_path_t hop,
                               const size_t* src,
                               size_t Nfs,
                               vector<size_t>& dest) const
{
    dest.resize(m_Lx*m_Ly);
    memcpy(&dest[0],src,m_Lx*m_Ly*sizeof(size_t));
    for(size_t h=0;h<hop.size();++h){
        dest[hop[h].second]=dest[hop[h].first];
        dest[hop[h].first]=Nfs;
    }
}

hop_path_t WaveFunction::get_hop_path(
                const size_t* src,
                const size_t* dest,
                size_t Nfs) const
{
    hop_path_t out;
    // find differences between the two states
    vector<size_t> part,hole;
    for(size_t s=0;s<m_Lx*m_Ly;++s){
        if(src[s]==Nfs && dest[s]!=Nfs){
            part.push_back(s);
        }
        if(src[s]!=Nfs && dest[s]==Nfs){
            hole.push_back(s);
        }
    }
    for(size_t h=0;h<part.size();++h)
        out.push_back(pair<size_t,size_t>(hole[h],part[h]));
    return out;
}

void WaveFunction::GetHops(vector<hop_path_t>& hopup,
                           vector<hop_path_t>& hopdo) const
{
    hopup.resize(m_excup.size());
    hopdo.resize(hopup.size());
    for(size_t s=0;s<hopup.size();++s){
        hopup[s]=hop_path_t();
        for(size_t h=0;h<m_excup[m_c_exc][s].size();++h){
            hopup[s].push_back(
                    pair<size_t,size_t>(
                        m_fockup[m_excup[m_c_exc][s][h].first],
                        m_excup[m_c_exc][s][h].second));
        }
        hopdo[s]=hop_path_t();
        for(size_t h=0;h<m_excdo[m_c_exc][s].size();++h){
            hopdo[s].push_back(
                    pair<size_t,size_t>(
                        m_fockdo[m_excdo[m_c_exc][s][h].first],
                        m_excdo[m_c_exc][s][h].second));
        }
    }
}

void WaveFunction::GetHop(size_t k,
                          hop_path_t& hopup,
                          hop_path_t& hopdo) const
{
    hopup.clear();
    hopdo.clear();
    for(size_t h=0;h<m_excup[m_c_exc][k].size();++h)
        hopup.push_back(pair<size_t,size_t>(
                    m_fockup[m_excup[m_c_exc][k][h].first],
                    m_excup[m_c_exc][k][h].second));
    for(size_t h=0;h<m_excdo[m_c_exc][k].size();++h)
        hopdo.push_back(pair<size_t,size_t>(
                    m_fockdo[m_excdo[m_c_exc][k][h].first],
                    m_excdo[m_c_exc][k][h].second));
}

void WaveFunction::hop(size_t khop)
{
    hop_path_t hopup,hopdo;
    GetHop(khop,hopup,hopdo);
    m_sign*=hop_sign(hopup,hopdo);
    // update m_fockup and m_fockdo
    for(size_t h=0;h<hopup.size();++h){
        m_fockup[m_fsup[hopup[h].first]]=m_Nfsup;
        m_fockup[hopup[h].second]=hopup[h].first;
        m_fsup[hopup[h].first]=hopup[h].second;
    }
    for(size_t h=0;h<hopdo.size();++h){
        m_fockdo[m_fsdo[hopdo[h].first]]=m_Nfsdo;
        m_fockdo[hopdo[h].second]=hopdo[h].first;
        m_fsdo[hopdo[h].first]=hopdo[h].second;
    }
    m_c_exc=khop;
}

string WaveFunction::Fock() const
{
    ostringstream out;
    out<<"up  : |";
    for(size_t f=0;f<m_Lx*m_Ly;++f){
        if(m_fockup[f]==m_Nfsup)
            out<<std::setw(2)<<"."<<" ";
        else
            out<<std::setw(2)<<m_fockup[f]<<" ";
    }
    out<<">"<<endl<<"down: |";
    for(size_t f=0;f<m_Lx*m_Ly;++f){
        if(m_fockdo[f]==m_Nfsdo)
            out<<std::setw(2)<<"."<<" ";
        else
            out<<std::setw(2)<<m_fockdo[f]<<" ";
    }
    out<<">"<<endl;
    return out.str();
}

void WaveFunction::save(FileManager* fm)
{
    int rank=0, size=1;
#ifdef USEMPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif//USEMPI
    if(size==1 || rank==1){ //not root process except if non-mpi
        hid_t wfile=fm->WriteSimple("WaveFunction");
        hsize_t dims[2]={m_fock_states_up.size(),m_Lx*m_Ly};
        int* dup=new int[m_fock_states_up.size()*m_Lx*m_Ly];
        int* ddo=new int[m_fock_states_up.size()*m_Lx*m_Ly];
        for(size_t s=0;s<m_fock_states_up.size();++s){
            memcpy(&dup[s*m_Lx*m_Ly],&m_fock_states_up[s][0],m_Lx*m_Ly*sizeof(int));
            memcpy(&ddo[s*m_Lx*m_Ly],&m_fock_states_do[s][0],m_Lx*m_Ly*sizeof(int));
        }
        H5LTmake_dataset_int(wfile,"/states_up",2,dims,dup);
        H5LTmake_dataset_int(wfile,"/states_do",2,dims,ddo);
        delete [] dup;
        delete [] ddo;
        H5Fclose(wfile);
    }
}

std::ostream & operator<<(std::ostream& out, const WaveFunction & wav)
{
    out<<"********** WaveFunction: *************"<<std::endl;
    out<<"global sign="<<wav.GetSign()<<std::endl;
    out<<"up spins by indices"<<std::endl;
    for(size_t k=0;k<wav.m_Nfsup;++k)
        out<<std::setw(2)<<k<<" ";
    out<<std::endl;
    for(size_t k=0;k<wav.m_Nfsup;++k)
        out<<"---";
    out<<std::endl;
    for(size_t k=0;k<wav.m_Nfsup;++k)
        out<<std::setw(2)<<wav.m_fsup[k]<<" ";
    out<<std::endl;
    out<<std::endl<<" Fock state for spins up:"<<std::endl;
    for(size_t f=0;f<wav.m_Lx*wav.m_Ly;++f){
        if(wav.m_fockup[f]==wav.m_Nfsup)
            out<<std::setw(2)<<"."<<" ";
        else
            out<<std::setw(2)<<wav.m_fockup[f]<<" ";
    }
    out<<">"<<std::endl;
    out<<"down spins by indices"<<std::endl;
    for(size_t k=0;k<wav.m_Nfsdo;++k)
        out<<std::setw(2)<<k<<" ";
    out<<std::endl;
    for(size_t k=0;k<wav.m_Nfsdo;++k)
        out<<"---";
    out<<std::endl;
    for(size_t k=0;k<wav.m_Nfsdo;++k)
        out<<std::setw(2)<<wav.m_fsdo[k]<<" ";
    out<<std::endl;
    out<<std::endl<<" Fock state for spins down:"<<std::endl;
    for(size_t f=0;f<wav.m_Lx*wav.m_Ly;++f){
        if(wav.m_fockdo[f]==wav.m_Nfsdo)
            out<<std::setw(2)<<"."<<" ";
        else
            out<<std::setw(2)<<wav.m_fockdo[f]<<" ";
    }
    out<<">"<<std::endl;
    out<<std::endl;
    return out;
}    
