#include "Jastrow.h"
#include "LatticeState.h"
#include "JastrowPotential.h"
#include "Lattice.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <cassert>

using namespace std;

Jastrow::Jastrow(const LatticeState* st,const JastrowPotential* pot)
    :Amplitude(st),
     m_jassum(0),
     m_jaspot(pot)
{
    if(pot){
        m_jassum_grad=vector<double>(pot->GetNParams(),0);
        m_jassum_hess=vector<double>(pow(pot->GetNParams(),2),0);
    }
}

void Jastrow::Init()
{
    m_jassum=0;
    m_jassum_grad.assign(m_jassum_grad.size(),0);
    m_jassum_hess.assign(m_jassum_hess.size(),0);
    vector<double> grad(m_jaspot->GetNParams(),0);
    vector<double> hess(pow(m_jaspot->GetNParams(),2),0);
    for(size_t fli=0;fli<m_latstate->GetNfl();++fli){
        for(size_t pi=0;pi<m_latstate->GetNpt()[fli];++pi){
            for(size_t flj=0;flj<m_latstate->GetNfl();++flj){
                for(size_t pj=0;pj<m_latstate->GetNpt()[flj];++pj){
                    uint_vec_t statei(4),statej(4);
                    m_latstate->Fock2QN(m_latstate->Getpt()[fli][pi],fli,statei);
                    m_latstate->Fock2QN(m_latstate->Getpt()[flj][pj],flj,statej);
                    m_jassum+=m_jaspot->Pot(statei,statej);
                    m_jaspot->PotGrad(statei,statej,grad);
                    m_jaspot->PotHess(statei,statej,hess);
                    for(size_t pa=0;pa<m_jaspot->GetNParams();++pa){
                        m_jassum_grad[pa]+=grad[pa];
                        for(size_t pb=0;pb<m_jaspot->GetNParams();++pb)
                            m_jassum_hess[pa*m_jaspot->GetNParams()+pb]+=
                                hess[pa*m_jaspot->GetNParams()+pb];
                    }
                }
            }
        }
    }
} 

size_t Jastrow::GetNParams() const
{
    return m_jaspot->GetNParams();
}

double Jastrow::Jas() const
{
    return exp(m_jassum);
}

void Jastrow::JasGrad(vector<double>& grad) const
{
    grad.resize(m_jaspot->GetNParams(),0);
    for(size_t p=0;p<m_jaspot->GetNParams();++p)
        grad[p]=exp(m_jassum)*m_jassum_grad[p];
}

void Jastrow::JasHess(vector<double>& hess) const
{
    hess.resize(pow(m_jaspot->GetNParams(),2),0);
    for(size_t pa=0;pa<m_jaspot->GetNParams();++pa)
        for(size_t pb=0;pb<m_jaspot->GetNParams();++pb)
            hess[pa*m_jaspot->GetNParams()+pb]=exp(m_jassum)*(
                                               m_jassum_grad[pa]*m_jassum_grad[pb]+
                                               m_jassum_hess[pa*m_jaspot->GetNParams()+pb]);
}

void Jastrow::VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                         std::vector<double>& js) const
{
    vector<double> dummy_grad,dummy_hess;
    virtual_update(rhop,js,dummy_grad,dummy_hess,false);
}

void Jastrow::VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                         std::vector<double>& js,
                         std::vector<double>& js_grad,
                         std::vector<double>& js_hess) const
{
    virtual_update(rhop,js,js_grad,js_hess,true);
}

void Jastrow::virtual_update(const std::vector<std::vector<hop_path_t> >& rhop,
                             std::vector<double>& js,
                             std::vector<double>& js_grad,
                             std::vector<double>& js_hess,
                             bool grad_hess) const
{
    size_t np=m_jaspot->GetNParams();
    vector<double> gradnewij(np);
    vector<double> gradoldij(np);
    vector<double> hessnewij(pow(np,2));
    vector<double> hessoldij(pow(np,2));
    vector<double> gradnewji(np);
    vector<double> gradoldji(np);
    vector<double> hessnewji(pow(np,2));
    vector<double> hessoldji(pow(np,2));
    js.assign(rhop.size(),m_jassum);
    js_grad.resize(rhop.size()*np);
    js_hess.resize(rhop.size()*np*np);
    for(size_t r=0;r<rhop.size();++r){
        memcpy(&js_grad[r*np],&m_jassum_grad[0],np*sizeof(size_t));
        memcpy(&js_hess[r*np*np],&m_jassum_hess[0],np*np*sizeof(size_t));
    }
    for(size_t hidx=0;hidx<rhop.size();++hidx){
        //list of particles which have moved.
        vector<uint_vec_t> moving_part(m_latstate->GetNfl());
        for(size_t flk=0;flk<m_latstate->GetNfl();++flk){
            for(size_t h=0;h<rhop[hidx][flk].size();++h){
                moving_part[flk].push_back(m_latstate->Getfs()[flk][rhop[hidx][flk][h].first]);
            }
        }
        for(size_t flk=0;flk<m_latstate->GetNfl();++flk){
            for(size_t k=0;k<rhop[hidx][flk].size();++k){
                uint_vec_t statekold(3),stateknew(3);
                m_latstate->Fock2QN(rhop[hidx][flk][k].first,flk,statekold);
                m_latstate->Fock2QN(rhop[hidx][flk][k].second,flk,stateknew);
                assert(statekold[0]==m_latstate->GetLattice()->GetVertices()[statekold[0]]->idx);
                for(size_t fli=0;fli<m_latstate->GetNfl();++fli){
                    for(size_t pti=0;pti<m_latstate->GetNpt()[fli];++pti){
                        if(find(moving_part[fli].begin(),moving_part[fli].end(),pti)==moving_part[fli].end()){
                            uint_vec_t statei(3);
                            m_latstate->Fock2QN(m_latstate->Getpt()[fli][pti],fli,statei);
                            js[hidx]+=m_jaspot->Pot(statei,stateknew)+
                                m_jaspot->Pot(stateknew,statei)-
                                m_jaspot->Pot(statei,statekold)-
                                m_jaspot->Pot(statekold,statei);
                            if(grad_hess){
                                m_jaspot->PotGrad(statei,stateknew,gradnewij);
                                m_jaspot->PotGrad(stateknew,statei,gradnewji);
                                m_jaspot->PotGrad(statei,statekold,gradoldij);
                                m_jaspot->PotGrad(statekold,statei,gradoldji);
                                m_jaspot->PotHess(statei,stateknew,hessnewij);
                                m_jaspot->PotHess(stateknew,statei,hessnewji);
                                m_jaspot->PotHess(statei,statekold,hessoldij);
                                m_jaspot->PotHess(statekold,statei,hessoldji);
                                for(size_t pa=0;pa<np;++pa){
                                    js_grad[hidx*np+pa]+=gradnewij[pa]+gradnewji[pa]-
                                                         gradoldij[pa]-gradoldji[pa];
                                    for(size_t pb=0;pb<np;++pb)
                                        js_hess[(hidx*np+pa)*np+pb]+=hessnewij[pa*np+pb]
                                                                    +hessnewji[pa*np+pb]
                                                                    -hessoldij[pa*np+pb]
                                                                    -hessoldji[pa*np+pb];
                                }
                            }
                        }
                    }
                }
                for(size_t fll=0;fll<m_latstate->GetNfl();++fll){
                    for(size_t l=0;l<rhop[hidx][fll].size();++l){
                        uint_vec_t statelold(4),statelnew(3);
                        m_latstate->Fock2QN(rhop[hidx][fll][l].first,fll,statelold);
                        m_latstate->Fock2QN(rhop[hidx][fll][l].second,fll,statelnew);
                        js[hidx]+=m_jaspot->Pot(stateknew,statelnew)-
                            m_jaspot->Pot(statekold,statelold);
                        if(grad_hess){
                            m_jaspot->PotGrad(stateknew,statelnew,gradnewij);
                            m_jaspot->PotGrad(statekold,statelold,gradoldij);
                            m_jaspot->PotHess(stateknew,statelnew,hessnewij);
                            m_jaspot->PotHess(statekold,statelold,hessoldij);
                            for(size_t pa=0;pa<np;++pa){
                                js_grad[hidx*np+pa]+=gradnewij[pa]-gradoldij[pa];
                                for(size_t pb=0;pb<np;++pb)
                                    js_hess[(hidx*np+pa)*np+pb]+=hessnewij[pa*np+pb]
                                                                -hessoldij[pa*np+pb];
                            }
                        }
                    }
                }
            }
        }
    }
    for(size_t j=0;j<js.size();++j){
        js[j]=exp(js[j]);
        for(size_t pa=0;pa<np;++pa){
            for(size_t pb=0;pb<np;++pb){
                js_hess[(j*np+pa)*np+pb]=js[j]*(js_grad[j*np+pa]*js_grad[j*np+pb]+
                                                js_hess[(j*np+pa)*np+pb]);
            }
        }
        for(size_t p=0;p<np;++p){
            js_grad[j*np+p]=js[j]*js_grad[j*np+p];
        }
    }
}

void Jastrow::Update(const vector<hop_path_t>& rhop)
{
    size_t np=m_jaspot->GetNParams();
    //list of particles which have moved.
    vector<uint_vec_t> moving_part(m_latstate->GetNfl());
    for(size_t flk=0;flk<m_latstate->GetNfl();++flk){
        for(size_t h=0;h<rhop[flk].size();++h){
            moving_part[flk].push_back(m_latstate->Getfs()[flk][rhop[flk][h].second]);
        }
    }
    vector<double> gradnewij(np);
    vector<double> gradoldij(np);
    vector<double> hessnewij(pow(np,2));
    vector<double> hessoldij(pow(np,2));
    vector<double> gradnewji(np);
    vector<double> gradoldji(np);
    vector<double> hessnewji(pow(np,2));
    vector<double> hessoldji(pow(np,2));
    for(size_t flk=0;flk<m_latstate->GetNfl();++flk){
        for(size_t k=0;k<rhop[flk].size();++k){
            uint_vec_t statekold(3),stateknew(3);
            m_latstate->Fock2QN(rhop[flk][k].first,flk,statekold);
            m_latstate->Fock2QN(rhop[flk][k].second,flk,stateknew);
            for(size_t fli=0;fli<m_latstate->GetNfl();++fli){
                for(size_t pti=0;pti<m_latstate->GetNpt()[fli];++pti){
                    //if fli.pti is not one of the moving particles:
                    if(find(moving_part[fli].begin(),moving_part[fli].end(),pti)==moving_part[fli].end()){
                        uint_vec_t statei(3);
                        m_latstate->Fock2QN(m_latstate->Getpt()[fli][pti],fli,statei);
                        m_jassum+=m_jaspot->Pot(statei,stateknew)+
                            m_jaspot->Pot(stateknew,statei)-
                            m_jaspot->Pot(statei,statekold)-
                            m_jaspot->Pot(statekold,statei);
                        m_jaspot->PotGrad(statei,stateknew,gradnewij);
                        m_jaspot->PotGrad(stateknew,statei,gradnewji);
                        m_jaspot->PotGrad(statei,statekold,gradoldij);
                        m_jaspot->PotGrad(statekold,statei,gradoldji);
                        m_jaspot->PotHess(statei,stateknew,hessnewij);
                        m_jaspot->PotHess(stateknew,statei,hessnewji);
                        m_jaspot->PotHess(statei,statekold,hessoldij);
                        m_jaspot->PotHess(statekold,statei,hessoldji);
                        for(size_t pa=0;pa<np;++pa){
                            m_jassum_grad[pa]+=gradnewij[pa]+gradnewji[pa]-
                                               gradoldij[pa]-gradoldji[pa];
                            for(size_t pb=0;pb<np;++pb)
                                m_jassum_hess[pa*np+pb]+=hessnewij[pa*np+pb]
                                                        +hessnewji[pa*np+pb]
                                                        -hessoldij[pa*np+pb]
                                                        -hessoldji[pa*np+pb];
                        }
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
                    m_jaspot->PotGrad(stateknew,statelnew,gradnewij);
                    m_jaspot->PotGrad(statekold,statelold,gradoldij);
                    m_jaspot->PotHess(stateknew,statelnew,hessnewij);
                    m_jaspot->PotHess(statekold,statelold,hessoldij);
                    for(size_t pa=0;pa<np;++pa){
                        m_jassum_grad[pa]+=gradnewij[pa]-gradoldij[pa];
                        for(size_t pb=0;pb<np;++pb)
                            m_jassum_hess[pa*np+pb]+=hessnewij[pa*np+pb]-hessoldij[pa*np+pb];
                    }
                }
            }
        }
    }
}
