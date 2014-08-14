#ifdef USEMPI
#include <mpi.h>
#endif
#include <vector>
#include "LatticeState.h"
#include "SquareLattice.h"
#include "StagFluxTransExciton.h"
#include "StagFluxLongExciton.h"
#include "StagFluxGroundState.h"
#include "SFpNpHxWaveFunction.h"
#include "SFpNpHxGroundState.h"
#include "SFpNxpHzWaveFunction.h"
#include "SFpNxpHzGroundState.h"
#include "SFpNxpHzExciton.h"
#include "SFpNpHxExciton.h"
#include "LatticeStepper.h"
#include "SlaterDeterminant.h"
#include "IdentityJastrow.h"
#include "Jastrow.h"
#include "NeelJastrowPotential.h"
#include "StagJastrowPotential.h"
#include "LogJastrowPotential.h"
#include "FileManager.h"
#include "ArgParse.h"
#include "RanGen.h"
#include <hdf5.h>
#include <hdf5_hl.h>

int main(int argc,char* argv[])
{
#ifdef USEMPI
    MPI_Init(&argc,&argv);
#endif
    int world_rank=0,world_size=1;
#ifdef USEMPI
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
#endif
    map<string,bool> bomap;
    map<string,size_t> simap;
    map<string,int> inmap;
    map<string,double> domap;
    map<string,string> stmap;
    simap["L"]=8;
    simap["steps"]=1000;
    simap["qx"]=0;
    simap["qy"]=0;
    inmap["seed"]=time(NULL);
    domap["phi"]=0.085;
    domap["neel"]=0.055;
    domap["mag_field"]=0.0;
    domap["jastrow"]=0.0;
    domap["phase_shift_x"]=1.0;
    domap["phase_shift_y"]=1.0;
    bomap["stagflux_wav"]=false;
    bomap["sfpnzphx_wav"]=false;
    bomap["sfpnxphz_wav"]=false;
    bomap["neel_jastrow"]=false;
    bomap["stag_jastrow"]=false;
    bomap["log_jastrow"]=false;
    stmap["spinstate"]="";
    stmap["channel"]="groundstate";
    int help=ArgParse(argc,argv,domap,inmap,simap,bomap,stmap);
    if(help){
#ifdef USEMPI
        MPI_Finalize();
#endif
        exit(0);
    }
    RanGen::srand(inmap["seed"]+world_rank);
    FileManager * fm=new FileManager(".","test_random_walk");
    domap["phi"]*=M_PI;
    for(map<string,bool>::iterator it=bomap.begin();it!=bomap.end();++it)
        fm->FileAttribute(it->first,it->second);
    for(map<string,size_t>::iterator it=simap.begin();it!=simap.end();++it)
        fm->FileAttribute(it->first,it->second);
    for(map<string,int>::iterator it=inmap.begin();it!=inmap.end();++it)
        fm->FileAttribute(it->first,it->second);
    for(map<string,double>::iterator it=domap.begin();it!=domap.end();++it)
        fm->FileAttribute(it->first,it->second);
    for(map<string,string>::iterator it=stmap.begin();it!=stmap.end();++it)
        fm->FileAttribute(it->first,it->second);
    fm->FileAttribute("gitversion",GIT_SHA1);
    if(world_rank==1 || world_size==1){
#ifdef USEMPI
        int calc_rank=0,calc_size=1;
        MPI_Comm_rank(fm->GetCalcComm(),&calc_rank);
        MPI_Comm_size(fm->GetCalcComm(),&calc_size);
#endif
        double rej=0;
        size_t L=simap["L"];
        vector<size_t> Q(2);
        double neel=domap["neel"];
        double phi=domap["phi"];
        double field=domap["mag_field"];
        vector<double> phase_shift(2);
        phase_shift[0]=domap["phase_shift_x"];
        phase_shift[1]=domap["phase_shift_y"];
        Q[0]=simap["qx"];
        Q[1]=simap["qy"];
        SquareLattice slat(L,L);
        LatticeState* sp(0);
        if(bomap["stagflux_wav"]){
            if(stmap["channel"]=="groundstate" || stmap["channel"]=="long"){
                sp=new LatticeState(fm,&slat,{L*L/2,L*L/2},{1,1});
            } else if(stmap["channel"]=="trans"){
                sp=new LatticeState(fm,&slat,{L*L/2+1,L*L/2-1},{1,1});
            }
        } else {
            sp=new LatticeState(fm,&slat,{L*L},{2});
        }
        WaveFunction* wav(0);
        if(bomap["stagflux_wav"]){
            if(stmap["channel"]=="groundstate"){
                wav=new StagFluxGroundState(fm,L,L,phi,neel,phase_shift);
            } else if(stmap["channel"]=="trans"){
                wav=new StagFluxTransExciton(fm,L,L,phi,neel,phase_shift,Q);
            } else if(stmap["channel"]=="long"){
                wav=new StagFluxLongExciton(fm,L,L,phi,neel,phase_shift,Q);
            }
        } else if(bomap["sfpnzphx_wav"]){
            if(stmap["channel"]=="groundstate"){
                wav=new SFpNpHxGroundState(fm,L,L,phi,neel,domap["hx"],phase_shift);
            } else {//trans and long are mixed
                wav=new SFpNpHxExciton(fm,L,L,phi,neel,domap["hx"],phase_shift,Q);
            }
        } else {
            if(stmap["channel"]=="groundstate"){
                wav=new SFpNxpHzGroundState(fm,L,L,phi,neel,domap["hz"],phase_shift);
            } else {//trans and long are mixed
                wav=new SFpNxpHzExciton(fm,L,L,phi,neel,domap["hz"],phase_shift,Q);
            }
        }
        SlaterDeterminant amp(sp,wav);
        Jastrow* jas=0;
        JastrowPotential* jaspot=0;
        if(bomap["neel_jastrow"]){
            jaspot=new NeelJastrowPotential(&slat,domap["jastrow"]);
            jas=new Jastrow(sp,jaspot);
        } else if(bomap["stag_jastrow"]){
            jaspot=new StagJastrowPotential(&slat,domap["jastrow"]);
            jas=new Jastrow(sp,jaspot);
        } else if(bomap["log_jastrow"]){
            jaspot=new LogJastrowPotential(&slat,domap["jastrow"]);
            jas=new Jastrow(sp,jaspot);
        } else {
            jas=new IdentityJastrow;
        }
        if(stmap["spinstate"]!=""){
            int calc_rank=0, calc_size=1;
#ifdef USEMPI
            MPI_Comm_rank(fm->GetCalcComm(),&calc_rank);
            MPI_Comm_size(fm->GetCalcComm(),&calc_size);
#endif
            vector<vector<int> > fockstates(sp->GetNfl());
            if(calc_rank==0){
                hid_t ifile=H5Fopen(stmap["spinstate"].c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
#ifdef USEMPI
                for(size_t fl=0;fl<sp->GetNfl();++fl){
                    fockstates[fl]=vector<int>(sp->GetNfs()[fl],0);
                    vector<int> rstate(sp->GetNfs()[fl]*calc_size);
                    for(size_t r=0;r<calc_size;++r){
                        H5LTread_dataset_int(ifile,
                                (string("/rank-")+to_string(r+1)+
                                 string("/flavour-")+to_string(fl)).c_str(),
                                &rstate[r*sp->GetNfs()[fl]]);
                    }
                    MPI_Scatter(&rstate[0],sp->GetNfs()[fl],MPI_INT,
                            &fockstates[fl][0],sp->GetNfs()[fl],MPI_INT,
                            0,fm->GetCalcComm());
                }
#else
                for(size_t fl=0;fl<sp->GetNfl();++fl){
                    fockstates[fl]=vector<int>(sp->GetNfs()[fl],0);
                    H5LTread_dataset_int(ifile,
                        (string("/rank-0/flavour-")+to_string(fl)).c_str(),
                        &fockstates[fl][0]);
                }
#endif
            } else {
#ifdef USEMPI
                for(size_t fl=0;fl<sp->GetNfl();++fl){
                    fockstates[fl]=vector<int>(sp->GetNfs()[fl],0);
                    MPI_Scatter(NULL,0,MPI_INT,
                            &fockstates[fl][0],sp->GetNfs()[fl],MPI_INT,
                            0,fm->GetCalcComm());
                }
#endif
            }
            sp->FockInit(fockstates);
            amp.Init();
            jas->Init();
        } else {
            vector<vector<size_t> > pop;
            if(bomap["stagflux_wav"]){
                if(stmap["channel"]=="trans"){
                    pop.push_back(vector<size_t>(1,L*L/2+1));
                    pop.push_back(vector<size_t>(1,L*L/2-1));
                } else {
                    pop.push_back(vector<size_t>(1,L*L/2));
                    pop.push_back(vector<size_t>(1,L*L/2));
                }
            } else {
                pop.push_back({L*L/2,L*L/2});
            }
            sp->RanInit(pop);
            amp.Init();
            int tc=0,failcount=10;
            while(amp.Amp()==0.0 && tc<failcount){
                sp->RanInit(pop);
                amp.Init();
                jas->Init();
                tc++;
            }
            if(tc==failcount){
                cerr<<"could not find a non-singular starting configuration."<<endl;
                cerr<<"last real space state is:"<<endl<<*sp<<endl;
                cerr<<"wf state is:"<<endl<<*wav<<endl;
#ifdef USEMPI
                MPI_Finalize();
#endif
                exit(1);
            }
        }
        BigDouble maxq=norm(amp.Amp());
        size_t kmax=wav->GetCurrentStateIndex();
        for(size_t k=0;k<wav->GetNExc();++k){
            wav->Hop(k);
            amp.Init();
            if(norm(amp.Amp())>maxq){
                kmax=k;
                maxq=norm(amp.Amp());
            }
        }
        wav->Hop(kmax);
        amp.Init();
        LatticeStepper step(sp,wav,&amp,jas);
        vector<bool> flip(sp->GetNfl(),sp->GetNfl()==1);
        step.SetFlavorFlip(flip);
        for(size_t s=0;s<simap["steps"];++s){
            BigDouble weight=step.weight();
            cout<<"weight: "<<weight<<endl;
            BigDouble tryweight=step.trystep();
            cout<<"trial weight: "<<tryweight<<endl;
            double ratio=(double)(tryweight/weight);
            double r=RanGen::uniform();
            if(ratio>1.0 || r<ratio){
                step.step();
                cout<<(double)(step.weight()/tryweight)<<endl;
            }
        }
    }
}
