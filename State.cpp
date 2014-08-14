#include "State.h"
#include "FileManager.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#ifdef USEMPI
#include <mpi.h>
#endif
#include <iomanip>

using namespace std;

State::State(FileManager* fm)
    : m_filemanager(fm),m_sign(1)
{}

State::~State()
{}

void State::InitPart(const vector<uint_vec_t>& st)
{
    m_part=st;
    for(size_t fl=0;fl<m_Nfl;++fl){
        for(size_t p=0;p<m_Npt[fl];++p){
            m_fock[fl][m_part[fl][p]]=p;
        }
    }
}

void State::InitFock(const vector<uint_vec_t>& st)
{
    for(size_t fl=0;fl<m_Nfl;++fl){
        size_t pi(0);
        for(size_t s=0;s<m_Nfs[fl];++s){
            if(st[fl][s]<m_Npt[fl]){
                m_fock[fl][s]=pi;
                m_part[fl][pi]=s;
                ++pi;
            }
        }
    }
}

const uint_vec_t& State::GetNfs() const
{
    return m_Nfs;
}

const uint_vec_t& State::GetNpt() const
{
    return m_Npt;
}

size_t State::GetNfl() const
{
    return m_Nfl;
}

const vector<uint_vec_t>& State::Getfs() const
{
    return m_fock;
}

const vector<uint_vec_t>& State::Getpt() const
{
    return m_part;
}

int State::GetSign() const
{
    return m_sign;
}

void State::Hop(const vector<hop_path_t>& hops)
{
    m_sign*=HopSign(hops);
    for(size_t fl=0;fl<m_Nfl;++fl){
        for(size_t h=0;h<hops[fl].size();++h){
            size_t part=m_fock[fl][hops[fl][h].first];
            m_fock[fl][hops[fl][h].first]=m_Npt[fl];
            m_fock[fl][hops[fl][h].second]=part;
            m_part[fl][part]=hops[fl][h].second;
        }
    }
}

int State::HopSign(const vector<hop_path_t>& hops) const
{
    // calculate sign.
    bool odd=false;
    size_t mn,mx;
    for(size_t fl=0;fl<m_Nfl;++fl){
        if(hops[fl].size()){
            mn=min(hops[fl][0].first,hops[fl][0].second)+1;
            mx=max(hops[fl][0].first,hops[fl][0].second);
            for(size_t s=mn;s<mx;++s)
                odd= (odd != (m_fock[fl][s]!=m_Npt[fl]));
            if(hops[fl].size()>1){
                uint_vec_t sf(m_fock[fl]);
                size_t part=sf[hops[fl][0].first]; 
                sf[hops[fl][0].first]=m_Npt[fl];
                sf[hops[fl][0].second]=part;
                for(size_t h=1;h<hops[fl].size();++h){
                    mn=min(hops[fl][h].first,hops[fl][h].second)+1;
                    mx=max(hops[fl][h].first,hops[fl][h].second);
                    for(size_t s=mn;s<mx;++s)
                        odd = ( odd != (sf[s] != m_Npt[fl]));
                    part=sf[hops[fl][h].first];
                    sf[part]=m_Npt[fl];
                    sf[hops[fl][h].second]=part;
                }
            }
        }
    }
    return (2*int(!odd)-1);
}

void State::Save()
{
    int rank=0, size=1;
#ifdef USEMPI
    MPI_Comm comm=m_filemanager->GetCalcComm();
    MPI_Comm_size(comm,&size);
    MPI_Comm_rank(comm,&rank);
#endif
    vector<vector<int> > st;
    for(size_t fl=0;fl<m_Nfl;++fl){
        st.push_back(vector<int>(m_Nfs[fl],0));
        for(size_t f=0;f<m_Nfs[fl];++f)
            st[fl][f]=int(m_fock[fl][f]!=m_Npt[fl]);
    }
    if(rank==0){
        hid_t ofile=m_filemanager->WriteSimple("FockState");
        for(size_t fl=0;fl<m_Nfl;++fl){
            hsize_t dim=m_Nfs[fl];
            vector<int> buf(m_Nfs[fl]*(size));
#ifdef USEMPI
            MPI_Gather(&st[fl][0],m_Nfs[fl],MPI_INT,&buf[0],m_Nfs[fl],MPI_INT,0,comm);
#else
            buf=st[fl];
#endif
            for(size_t r=0;r<size;++r){
                herr_t status;
                hid_t grp;
                int absrankid=0;
#ifdef USEMPI
                absrankid=r+1;
#endif
                bool grpexists=H5Lexists(ofile,(string("/rank-")+to_string(absrankid)).c_str(),H5P_DEFAULT);
                if(grpexists){
                    grp=H5Gopen(ofile,(string("/rank-")+to_string(absrankid)).c_str(),H5P_DEFAULT);
                } else {
                    grp=H5Gcreate(ofile,(string("rank-")+to_string(absrankid)).c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
                }
                H5LTmake_dataset_int(grp,(string("flavour-")+to_string(fl)).c_str(),1,&dim,&buf[r*m_Nfs[fl]]);
                H5Gclose(grp);
            }
        }
        H5Fclose(ofile);
    } else {
#ifdef USEMPI
        for(size_t fl=0;fl<m_Nfl;++fl){
            MPI_Gather(&st[fl][0],m_Nfs[fl],MPI_INT,NULL,0,MPI_INT,0,comm);
        }
#endif
    }
}

std::ostream& operator<<(std::ostream& out,
                         const State& wav)
{
    out<<"********** State: *************"<<std::endl;
    out<<"global sign="<<wav.GetSign()<<std::endl;
    for(size_t fl=0;fl<wav.m_Nfl;++fl){
        out<<"flavour "<<fl<<" by particle indices"<<std::endl;
        for(size_t k=0;k<wav.m_Npt[fl];++k)
            out<<std::setw(2)<<k<<" ";
        out<<std::endl;
        for(size_t k=0;k<wav.m_Npt[fl];++k)
            out<<"---";
        out<<std::endl;
        for(size_t k=0;k<wav.m_Npt[fl];++k)
            out<<std::setw(2)<<wav.m_part[fl][k]<<" ";
        out<<std::endl;
        out<<std::endl<<"flavour "<<fl<<" in Fock space:"<<std::endl;
        for(size_t f=0;f<wav.m_Nfs[fl];++f){
            if(wav.m_fock[fl][f]==wav.m_Npt[fl])
                out<<std::setw(2)<<"."<<" ";
            else
                out<<std::setw(2)<<wav.m_fock[fl][f]<<" ";
        }
        out<<">"<<std::endl;
    }
    out<<std::endl;
    return out;
}

//################## protected functions ####################

void State::build_base(const uint_vec_t& Npt,
                         const uint_vec_t& Nfs)
{
    m_Nfl=Npt.size();
    m_part=vector<uint_vec_t>(m_Nfl);
    m_fock=vector<uint_vec_t>(m_Nfl);
    m_Nfs=Nfs;
    m_Npt=Npt;
    for(size_t fl=0;fl<m_Nfl;++fl){
        m_part[fl]=uint_vec_t(Npt[fl],Nfs[fl]);
        m_fock[fl]=uint_vec_t(Nfs[fl],Npt[fl]);
    }
}

void State::get_hop_path(const uint_vec_t& src,
                           const uint_vec_t& dest,
                           const size_t& Npt,
                           hop_path_t& hop) const
{
    // find differences between the two states
    hop.clear();
    uint_vec_t part,hole;
    for(size_t s=0;s<src.size();++s){
        if(src[s]==Npt && dest[s]!=Npt){
            part.push_back(s);
        }
        if(src[s]!=Npt && dest[s]==Npt){
            hole.push_back(s);
        }
    }
    for(size_t h=0;h<part.size();++h)
        hop.push_back(pair<size_t,size_t>(hole[h],part[h]));
}

void State::get_hop_state(const hop_path_t& hop,
                            const uint_vec_t& src,
                            const size_t& Npt,
                            uint_vec_t& dest) const
{
    dest=src;
    for(size_t h=0;h<hop.size();++h){
        // No check whether hop is allowed (Pauli)
        dest[hop[h].second]=dest[hop[h].first];
        dest[hop[h].first]=Npt;
    }
}
