#include "State_1.h"
#include <iomanip>

using namespace std;

State_1::State_1()
    : m_sign(1)
{}

State_1::~State_1()
{}

void State_1::Init(vector<uint_vec_t> st)
{
    m_part=st;
    for(size_t fl=0;fl<m_Nfl;++fl){
        for(size_t p=0;p<m_Npt[fl];++p){
            m_fock[fl][m_part[fl][p]]=p;
        }
    }
}

const uint_vec_t& State_1::GetNfs() const
{
    return m_Nfs;
}

const uint_vec_t& State_1::GetNpt() const
{
    return m_Npt;
}

const vector<uint_vec_t>& Getfs() const
{
    return m_fock;
}

const vector<uint_vec_t>& Getpt() const
{
    return m_part;
}

int State_1::GetSign() const
{
    return m_sign;
}

void State_1::Hop(const vector<hop_path_t>& hops)
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

int State_1::HopSign(const vector<hop_path_t>& hops) const
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

std::ostream& operator<<(std::ostream& out,
                         const State_1& wav)
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

void State_1::build_base(const uint_vec_t& Npt,
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

void State_1::get_hop_path(const uint_vec_t& src,
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

void State_1::get_hop_state(const hop_path_t& hop,
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
