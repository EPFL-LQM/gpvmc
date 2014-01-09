#include "WaveFunction_1.h"
#include "FileManager.h"
#include <hdf5_hl.h>
#include <hdf5.h>
#include <cstring>
#ifdef USEMPI
#include <mpi.h>
#endif

using namespace std;

WaveFunction_1::WaveFunction_1()
    : m_sign(1)
{}

WaveFunction_1::WaveFunction_1(std::vector<size_t> Nby,
                               std::vector<size_t> Nfs)
    : m_sign(1)
{
    // Lx and Ly specify real space lattice
    // size which generate a complete
    // basis. The Fock space size is thus 2*Lx*Ly.
    build_base(Nby,Nfs);
}

void WaveFunction_1::build_base(vector<size_t> Nby,vector<size_t> Nfs)
{
    m_fs=vector<size_t*>(Nby.size());
    m_fock=vector<size_t*>(Nby.size());
    m_cache=vector<complex<double>* >(Nby.size());
    m_exc=vector<vector<vector<hop_path_t> > >(Nby.size());
    m_fock_states=vector<vector<vector<int> > >(Nby.size());
    m_Nflav=Nby.size();
    m_Nfs=Nfs;
    m_Nby=Nby;
    for(size_t flav=0;flav<m_Nflav;++flav){
        m_fs[flav]=new size_t[Nby[flav]];
        m_fock[flav]=new size_t[m_Nfs[flav]];
        m_cache[flav]=new complex<double>[int(std::pow(m_Nfs[flav],2))];
        for(size_t s=0;s<m_Nfs[flav];++s){
            m_fock[flav][s]=m_Nby[flav];
        }
    }
}

WaveFunction_1::~WaveFunction_1()
{
    for(size_t flav=0;flav<m_Nflav;++flav){
        delete [] m_fs[flav];
        delete [] m_fock[flav];
        delete [] m_cache[flav];
    }
}

void WaveFunction_1::init_matrices()
{
    for(size_t flav=0;flav<m_Nflav;++flav){
        for(size_t f=0;f<m_Nfs[flav];++f){
            for(size_t R=0;R<m_Nfs[flav];++R){
                m_cache[flav][f*m_Nfs[flav]+R]=
                    this->matrix_element(f,R,flav);
            }
        }
    }
}

int WaveFunction_1::HopSign(const vector<hop_path_t>& hop) const
{
    // calculate sign.
    bool odd=false;
    size_t mn,mx;
    for(size_t flav=0;flav<m_Nflav;++flav){
        if(hop[flav].size()){
            mn=min(hop[flav][0].second,m_fs[flav][hop[flav][0].first])+1;
            mx=max(hop[flav][0].second,m_fs[flav][hop[flav][0].first]);
            for(size_t s=mn;s<mx;++s)
                odd= (odd != (m_fock[flav][s]!=m_Nby[flav]));
            if(hop[flav].size()>1){
                size_t* sf=new size_t[m_Nfs[flav]];
                std::memcpy(sf,m_fock[flav],m_Nfs[flav]*sizeof(size_t));
                sf[m_fs[flav][hop[flav][0].first]]=m_Nby[flav];
                sf[hop[flav][0].second]=hop[flav][0].first;
                for(size_t h=1;h<hop[flav].size();++h){
                    mn=min(hop[flav][h].second,m_fs[flav][hop[flav][h].first])+1;
                    mx=max(hop[flav][h].second,m_fs[flav][hop[flav][h].first]);
                    for(size_t s=mn;s<mx;++s)
                        odd = ( odd != (sf[s] != m_Nby[flav]));
                    sf[m_fs[flav][hop[flav][h].first]]=m_Nby[flav];
                    sf[hop[flav][h].second]=hop[flav][h].first;
                }
                delete [] sf;
            }
        }
    }
    return (2*int(!odd)-1);
}

void WaveFunction_1::AddState(const vector<const size_t*> f)
{
    if(!m_exc[0].size()){
        for(size_t flav=0;flav<m_Nflav;++flav){
            size_t sti=0;
            for(size_t s=0;s<m_Nfs[flav];++s){
                if(f[flav][s]!=m_Nby[flav]){
                    m_fock[flav][s]=sti;
                    m_fs[flav][sti]=s;
                    ++sti;
                }
            }
        }
        m_sign=1;
        m_c_exc=0;
    }
    hop_path_t em;// empty hop
    // span hop path matrix with empty hops
    for(size_t flav=0;flav<m_Nflav;++flav){
        for(size_t i=0;i<m_exc[flav].size();++i)
            m_exc[flav][i].push_back(em);
        m_exc[flav].push_back(
                vector<hop_path_t>(m_exc[flav].size()+1,em));
    }
    for(size_t flav=0;flav<m_Nflav;++flav){
        vector<size_t> state(m_Nfs[flav],0);
        for(size_t i=0;i<m_exc[flav].size()-1;++i){
            get_hop_state(m_exc[flav][0][i],m_fock[flav],m_Nby[flav],m_Nfs[flav],state);
            m_exc[flav].back()[i]=get_hop_path(
                    f[flav],&state[0],m_Nby[flav],m_Nfs[flav]);
            m_exc[flav][i].back()=get_hop_path(
                    &state[0],f[flav],m_Nby[flav],m_Nfs[flav]);
        }
        m_fock_states[flav].push_back(vector<int>(m_Nfs[flav]));
        for(size_t ff=0;ff<m_Nfs[flav];++ff){
            m_fock_states[flav].back()[ff]=(f[flav][ff]!=m_Nby[flav]);
        }
    }
}

void WaveFunction_1::get_hop_state(hop_path_t hop,
                                   const size_t* src,
                                   size_t Nby,
                                   size_t Nfs,
                                   vector<size_t>& dest) const
{
    dest.resize(Nfs);
    memcpy(&dest[0],src,Nfs*sizeof(size_t));
    for(size_t h=0;h<hop.size();++h){
        dest[hop[h].second]=dest[hop[h].first];
        dest[hop[h].first]=Nby;
    }
}

hop_path_t WaveFunction_1::get_hop_path(const size_t* src,
                                          const size_t* dest,
                                          size_t Nby,
                                          size_t Nfs) const
{
    hop_path_t out;
    // find differences between the two states
    vector<size_t> part,hole;
    for(size_t s=0;s<Nfs;++s){
        if(src[s]==Nby && dest[s]!=Nby){
            part.push_back(s);
        }
        if(src[s]!=Nby && dest[s]==Nby){
            hole.push_back(s);
        }
    }
    for(size_t h=0;h<part.size();++h)
        out.push_back(pair<size_t,size_t>(hole[h],part[h]));
    return out;
}

void WaveFunction_1::GetHops(std::vector<std::vector<hop_path_t> >& hop) const
{
    for(size_t flav=0;flav<m_Nflav;++flav){
        hop[flav].resize(m_exc[flav].size());
        for(size_t s=0;s<hop[flav].size();++s){
            for(size_t h=0;h<m_exc[flav][m_c_exc][s].size();++h){
                hop[flav][s].push_back(
                        pair<size_t,size_t>(
                            m_fock[flav][m_exc[flav][m_c_exc][s][h].first],
                            m_exc[flav][m_c_exc][s][h].second));
            }
        }
    }
}

void WaveFunction_1::GetHop(size_t k, std::vector<hop_path_t>& hop) const
{
    for(size_t flav=0;flav<m_Nflav;++flav){
        hop[flav].clear();
        for(size_t h=0;h<m_exc[flav][m_c_exc][k].size();++h){
            hop[flav].push_back(
                    pair<size_t,size_t>(
                        m_fock[flav][m_exc[flav][m_c_exc][k][h].first],
                        m_exc[flav][m_c_exc][k][h].second));
        }
    }
}
void WaveFunction_1::Hop(size_t khop)
{
    vector<hop_path_t> hop(m_Nflav);
    GetHop(khop,hop);
    m_sign*=HopSign(hop);
    for(size_t flav=0;flav<m_Nflav;++flav){
        for(size_t h=0;h<hop[flav].size();++h){
            m_fock[flav][m_fs[flav][hop[flav][h].first]]=m_Nby[flav];
            m_fock[flav][hop[flav][h].second]=hop[flav][h].first;
            m_fs[flav][hop[flav][h].first]=hop[flav][h].second;
        }
    }
    m_c_exc=khop;
}

string WaveFunction_1::Fock() const
{
    ostringstream out;
    for(size_t flav=0;flav<m_Nflav;++flav){
        out<<"flavor="<<flav<<"  : |";
        for(size_t f=0;f<m_Nfs[flav];++f){
            if(m_fock[flav][f]==m_Nby[flav])
                out<<std::setw(2)<<"."<<" ";
            else
                out<<std::setw(2)<<m_fock[flav][f]<<" ";
        }
        out<<">"<<endl;
    }
    return out.str();
}

void WaveFunction_1::Save(FileManager* fm)
{
    int rank=0, size=1;
#ifdef USEMPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif//USEMPI
    if(size==1 || rank==1){ //not root process except if non-mpi
        hid_t wfile=fm->WriteSimple("WaveFunction");
        for(size_t flav=0;flav<m_Nflav;++flav){
            hsize_t dims[2]={m_fock_states[flav].size(),m_Nfs[flav]};
            int* dd=new int[m_fock_states[flav].size()*m_Nfs[flav]];
            for(size_t s=0;s<m_fock_states[flav].size();++s)
                memcpy(&dd[s*m_Nfs[flav]],&m_fock_states[flav][s][0],m_Nfs[flav]*sizeof(int));
            ostringstream sn;
            sn<<"states_"<<flav;
            H5LTmake_dataset_int(wfile,sn.str().c_str(),2,dims,dd);
            delete [] dd;
        }
        H5Fclose(wfile);
    }
}

std::ostream & operator<<(std::ostream& out, const WaveFunction_1 & wav)
{
    out<<"********** WaveFunction: *************"<<std::endl;
    out<<"global sign="<<wav.GetSign()<<std::endl;
    for(size_t flav=0;flav<wav.m_Nflav;++flav){
        out<<"spin flavour "<<flav<<" by indices"<<std::endl;
        for(size_t k=0;k<wav.m_Nby[flav];++k)
            out<<std::setw(2)<<k<<" ";
        out<<std::endl;
        for(size_t k=0;k<wav.m_Nby[flav];++k)
            out<<"---";
        out<<std::endl;
        for(size_t k=0;k<wav.m_Nby[flav];++k)
            out<<std::setw(2)<<wav.m_fs[flav][k]<<" ";
        out<<std::endl;
        out<<std::endl<<" Fock state for spin flavour "<<flav<<":"<<std::endl;
        for(size_t f=0;f<wav.m_Nfs[flav];++f){
            if(wav.m_fock[flav][f]==wav.m_Nby[flav])
                out<<std::setw(2)<<"."<<" ";
            else
                out<<std::setw(2)<<wav.m_fock[flav][f]<<" ";
        }
        out<<">"<<std::endl;
    }
    out<<std::endl;
    return out;
}
