#include "Jastrow.h"
#include "LatticeState.h"
#include "JastrowPotential.h"
#include <cmath>

using namespace std;

Jastrow::Jastrow(const LatticeState* st,const JastrowPotential* pot)
    :m_jassum(0),m_latstate(st),m_jaspot(pot)
{}

void Jastrow::Init()
{
    m_jassum=0;
    for(size_t fli=0;fli<m_latstate->GetNfl();++fli){
        for(size_t pi=0;pi<m_latstate->GetNpt()[fli];++pi){
            for(size_t flj=0;flj<m_latstate->GetNfl();++flj){
                for(size_t pj=0;pj<m_latstate->GetNpt()[flj];++pj){
                    uint_vec_t statei(4),statej(4);
                    m_latstate->Fock2QN(m_latstate->Getpt()[fli][pi],fli,statei);
                    m_latstate->Fock2QN(m_latstate->Getpt()[flj][pj],flj,statej);
                    m_jassum+=m_jaspot->Pot(statei,statej);
                }
            }
        }
    }
} 

double Jastrow::Jas() const
{
    return exp(m_jassum);
}

void Jastrow::VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                         std::vector<double>& js) const
{
    js.assign(rhop.size(),m_jassum);
    for(size_t hidx=0;hidx<rhop.size();++hidx){
        for(size_t flk=0;flk<m_latstate->GetNfl();++flk){
            for(size_t k=0;k<rhop[hidx][flk].size();++k){
                uint_vec_t statekold(3),stateknew(3);
                m_latstate->Fock2QN(rhop[hidx][flk][k].first,flk,statekold);
                m_latstate->Fock2QN(rhop[hidx][flk][k].second,flk,stateknew);
                for(size_t fli=0;fli<m_latstate->GetNfl();++fli){
                    for(size_t pti=0;pti<m_latstate->GetNpt()[fli];++pti){
                        if(fli!=flk || pti!=m_latstate->Getfs()[flk][rhop[hidx][flk][k].first]){
                            uint_vec_t statei(3);
                            m_latstate->Fock2QN(m_latstate->Getpt()[fli][pti],fli,statei);
                            js[hidx]+=m_jaspot->Pot(statei,stateknew)+
                                m_jaspot->Pot(stateknew,statei)-
                                m_jaspot->Pot(statei,statekold)-
                                m_jaspot->Pot(statekold,statei);
                        }
                    }
                }
                for(size_t fll=0;fll<m_latstate->GetNfl();++fll){
                    for(size_t l=0;l<rhop[hidx][fll].size();++l){
                        uint_vec_t statelold(3),statelnew(3);
                        m_latstate->Fock2QN(rhop[hidx][fll][l].first,fll,statelold);
                        m_latstate->Fock2QN(rhop[hidx][fll][l].second,fll,statelnew);
                        js[hidx]+=m_jaspot->Pot(stateknew,statelnew)-
                            m_jaspot->Pot(statekold,statelold);
                    }
                }
            }
        }
    }
    for(size_t j=0;j<js.size();++j)
        js[j]=exp(js[j]);
}

void Jastrow::Update(const vector<hop_path_t>& rhop)
{
    for(size_t flk=0;flk<m_latstate->GetNfl();++flk){
        for(size_t k=0;k<rhop[flk].size();++k){
            uint_vec_t statekold(3),stateknew(3);
            m_latstate->Fock2QN(rhop[flk][k].first,flk,statekold);
            m_latstate->Fock2QN(rhop[flk][k].second,flk,stateknew);
            for(size_t fli=0;fli<m_latstate->GetNfl();++fli){
                for(size_t pti=0;pti<m_latstate->GetNpt()[fli];++pti){
                    if(fli!=flk || pti!=m_latstate->Getfs()[flk][rhop[flk][k].second]){
                        uint_vec_t statei(3);
                        m_latstate->Fock2QN(m_latstate->Getpt()[fli][pti],fli,statei);
                        m_jassum+=m_jaspot->Pot(statei,stateknew)+
                            m_jaspot->Pot(stateknew,statei)-
                            m_jaspot->Pot(statei,statekold)-
                            m_jaspot->Pot(statekold,statei);
                    }
                }
            }
            for(size_t fll=0;fll<m_latstate->GetNfl();++fll){
                for(size_t l=0;l<rhop[fll].size();++l){
                    uint_vec_t statelold(3),statelnew(3);
                    m_latstate->Fock2QN(rhop[fll][l].first,fll,statelold);
                    m_latstate->Fock2QN(rhop[fll][l].second,fll,statelnew);
                    m_jassum+=m_jaspot->Pot(stateknew,statelnew)-
                        m_jaspot->Pot(statekold,statelold);
                }
            }
        }
    }
}
