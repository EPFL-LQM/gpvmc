#include <iostream>
#include <map>
#include <ctime>
#include <cmath>
#include "FileManager.h"
#include "LatticeState.h"
#include "SquareLattice.h"
#include "Jastrow.h"
#include "IdentityJastrow.h"
#include "NeelJastrowPotential.h"
#include "StagJastrowPotential.h"
#include "LogJastrowPotential.h"
#include "ArgParse.h"
#include "RanGen.h"
#include "defs.h"
#ifdef USEMPI
#include <mpi.h>
#endif

using namespace std;

bool uint_vec_comp(const uint_vec_t& a, const uint_vec_t& b){
    if(a.size()!=b.size())
        return false;
    return equal(a.begin(),a.end(),b.begin());
}


int main(int argc,char* argv[])
{
#ifdef USEMPI
    MPI_Init(&argc,&argv);
#endif
    map<string,bool> bomap;
    map<string,size_t> simap;
    map<string,int> inmap;
    map<string,double> domap;
    map<string,string> stmap;
    simap["L"]=8;
    simap["N_hop_tests"]=10;
    simap["Max_hops"]=10;
    inmap["seed"]=time(NULL);
    bomap["Sztot_conserved"]=false;
    bomap["Sztot_non_zero_init"]=false;
    bomap["Neel_jastrow"]=false;
    bomap["Stag_jastrow"]=false;
    bomap["Log_jastrow"]=false;
    bomap["Neel_init"]=false;
    bomap["grad_hessian_test"]=false;
    domap["neel_field"]=1.0;
    int help=ArgParse(argc,argv,domap,inmap,simap,bomap,stmap);
    if(help){
#ifdef USEMPI
        MPI_Finalize();
#endif
        exit(0);
    }
    int rank=0;
#ifdef USEMPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
    FileManager fm(".","test_jastrow");
    RanGen::srand(inmap["seed"]+100*rank);
    size_t L=simap["L"];
    SquareLattice slat(L,L);
    size_t passed(0);
    size_t Nfl=1;
    size_t Nfs=2*L*L;
    size_t Npt=L*L;
    if(bomap["Sztot_conserved"]){
        Nfl=2;
        Npt=L*L/2;
        Nfs=L*L;
    }
    vector<vector<hop_path_t> > rhop(simap["N_hop_tests"],vector<hop_path_t>(Nfl));
    for(size_t tix=0;tix<simap["N_hop_tests"];++tix){
        LatticeState* st;
        if(bomap["Sztot_conserved"])
            st=new LatticeState(&fm,&slat,{L*L/2,L*L/2},{1,1});
        else
            st=new LatticeState(&fm,&slat,{L*L},{2});
        Jastrow* jas=0;
        JastrowPotential* jaspot=0;
        if(bomap["Neel_jastrow"]){
            jaspot=new NeelJastrowPotential(&slat,domap["neel_field"]);
            jas=new Jastrow(st,jaspot);
        } else if(bomap["Stag_jastrow"]){
            jaspot=new StagJastrowPotential(&slat,domap["neel_field"]);
            jas=new Jastrow(st,jaspot); 
        } else if(bomap["Log_jastrow"]){
            jaspot=new LogJastrowPotential(&slat,domap["neel_field"]);
            jas=new Jastrow(st,jaspot);
        } else {
            jas=new IdentityJastrow;
        }
        if(bomap["Neel_init"]){
            vector<uint_vec_t> fst(Nfl);
            for(size_t fl=0;fl<Nfl;++fl){
                fst[fl]=uint_vec_t(Nfs,Npt);
            }
            for(size_t v=0;v<slat.GetNv();++v){
                bool even=(slat.GetVertices()[v]->uc[0]+slat.GetVertices()[v]->uc[1])%2;
                if(even){
                    if(bomap["Sztot_conserved"])
                        fst[0][v]=0;
                    else
                        fst[0][v*2]=0;
                } else {
                    if(bomap["Sztot_conserved"])
                        fst[1][v]=0;
                    else
                        fst[0][v*2+1]=0;
                }
            }
            st->InitFock(fst);
        } else {
            vector<vector<size_t> > pop;
            if(bomap["Sztot_conserved"]){
                pop.push_back(vector<size_t>(1,L*L/2));
                pop.push_back(vector<size_t>(1,L*L/2));
            } else {
                if(bomap["Sztot_non_zero_init"]){
                    int Nup=(size_t)(RanGen::uniform()*L*L);
                    pop.push_back({size_t(Nup),size_t(int(L*L)-Nup)});
                } else {
                    pop.push_back({L*L/2,L*L/2});
                }
            }
            st->RanInit(pop);
        }
        jas->Init();
        cout<<"Initial state:"<<*st<<endl;
        cout<<"jastrow: "<<jas->Jas()<<endl;
        bool ok=false;
        size_t f1,f2;
        const Vertex* v1,*v2;
        while(!ok){
            size_t vx1=RanGen::uniform()*st->GetNsites();
            size_t vx2=RanGen::uniform()*st->GetNsites();
            v1=st->GetLattice()->GetVertices()[vx1];
            v2=st->GetLattice()->GetVertices()[vx2];
            if(vx1==vx2) continue;
            vector<uint_vec_t> s1,s2;
            st->GetLatOc(vx1,s1);
            st->GetLatOc(vx2,s2);
            if(equal(s1.begin(),s1.end(),s2.begin(),uint_vec_comp)) continue;
            size_t f1=max_element(s1.begin(),s1.end(),
                          [](const uint_vec_t& a, const uint_vec_t& b)
                          {return a.size()<b.size();})-s1.begin();
            size_t f2=max_element(s2.begin(),s2.end(),
                          [](const uint_vec_t& a, const uint_vec_t& b)
                          {return a.size()<b.size();})-s2.begin();
            rhop[tix][f1].push_back(hop_t(v1->idx*st->GetNifs()[f1]+s1[f1][0],
                                          v2->idx*st->GetNifs()[f2]+s1[f1][0]));
            rhop[tix][f2].push_back(hop_t(v2->idx*st->GetNifs()[f2]+s2[f2][0],
                                          v1->idx*st->GetNifs()[f1]+s2[f2][0]));
            ok=true;

        }
        vector<double> vjs(1);
        vector<double> vjs_grad(1*jas->GetNParams());
        vector<double> vjs_hess(1*pow(jas->GetNParams(),2));
        jas->VirtUpdate(vector<vector<hop_path_t> >(rhop.begin()+tix,rhop.begin()+tix+1),vjs,vjs_grad,vjs_hess);
        st->Hop(rhop[tix]);
        jas->Update(rhop[tix]);
        double ujas;
        vector<double> ujas_grad(1*jas->GetNParams());
        vector<double> ujas_hess(1*pow(jas->GetNParams(),2));
        vector<double> ijas_grad(1*jas->GetNParams());
        vector<double> ijas_hess(1*pow(jas->GetNParams(),2));
        ujas=jas->Jas();
        jas->JasGrad(ujas_grad);
        jas->JasHess(ujas_hess);
        jas->Init();
        jas->JasGrad(ijas_grad);
        jas->JasHess(ijas_hess);
        if(!bomap["grad_hessian_test"]){
            if(abs(ujas-jas->Jas())<1e-6 && abs(vjs[0]-jas->Jas())<1e-6){
                ++passed;
                cout<<"succeeded with hops:"<<endl;
                for(size_t fl=0;fl<Nfl;++fl){
                    cout<<"flavour "<<fl<<endl;
                    for(size_t h=0;h<rhop[tix][fl].size();++h){
                        cout<<"    ("<<rhop[tix][fl][h].first<<","<<rhop[tix][fl][h].second<<")"<<endl;
                    }
                }
                cout<<"VirtUpdate result: "<<vjs[0]<<" Update result: "<<ujas<<" Init result: "<<jas->Jas()<<endl;
            } else {
                cout<<"failed with hops:"<<endl;
                for(size_t fl=0;fl<Nfl;++fl){
                    cout<<"flavour "<<fl<<endl;
                    for(size_t h=0;h<rhop[tix][fl].size();++h){
                        cout<<"    ("<<rhop[tix][fl][h].first<<","<<rhop[tix][fl][h].second<<")"<<endl;
                    }
                }
                cout<<"VirtUpdate result: "<<vjs[0]<<" Update result: "<<ujas<<" Init result: "<<jas->Jas()<<endl;
                cout<<"Final state is"<<*st<<endl;
            }
        } else if(jaspot){
            double grad_vidiff=0,grad_uidiff=0,hess_vidiff=0,hess_uidiff=0;
            for(size_t pa=0;pa<jaspot->GetNParams();++pa){
                grad_vidiff+=abs(vjs_grad[pa]-ijas_grad[pa]);
                grad_uidiff+=abs(ujas_grad[pa]-ijas_grad[pa]);
                for(size_t pb=0;pb<jaspot->GetNParams();++pb){
                    hess_vidiff+=abs(vjs_hess[pa*jaspot->GetNParams()+pb]-ijas_hess[pa*jaspot->GetNParams()+pb]);
                    hess_uidiff+=abs(ujas_hess[pa*jaspot->GetNParams()+pb]-ijas_hess[pa*jaspot->GetNParams()+pb]);
                }
            }
            size_t np=jaspot->GetNParams();
            cout<<"Gradients are (virt, update, init):"<<endl;
            for(size_t pa=0;pa<jaspot->GetNParams();++pa)
                cout<<vjs_grad[pa]<<"\t|\t"<<ujas_grad[pa]<<"\t|\t"<<ijas_grad[pa]<<endl;
            cout<<"Hessians are (virt, update, init):"<<endl;
            for(size_t pa=0;pa<jaspot->GetNParams();++pa){
                for(size_t pb=0;pb<jaspot->GetNParams();++pb)
                    cout<<vjs_hess[pa*np+pb]<<" ";
                cout<<"\t|\t";
                for(size_t pb=0;pb<jaspot->GetNParams();++pb)
                    cout<<ujas_hess[pa*np+pb]<<" ";
                cout<<"\t|\t";
                for(size_t pb=0;pb<jaspot->GetNParams();++pb)
                    cout<<ijas_hess[pa*np+pb]<<" ";
                cout<<endl;
            }
            if(grad_vidiff<1e-6 && grad_uidiff<1e-6 && hess_vidiff<1e-6 && hess_uidiff<1e-6)
                ++passed;
        }
        delete jas;
        delete jaspot;
        delete st;
    }
    cout<<int(double(passed)/simap["N_hop_tests"]*100)<<"% passed"<<endl;
    if(passed == simap["N_hop_tests"])
        passed=0;
    else
        passed=1;
#ifdef USEMPI
    MPI_Finalize();
#endif
    return passed;
}

