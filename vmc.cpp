#ifdef USEMPI
#include <mpi.h>
#endif
#include <ctime>
#include <vector>
#include "MetroMC.h"
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
#include "ProjHeis.h"
#include "StatSpinStruct.h"
#include "StagMagn.h"
#include "StagMagnZ.h"
#include "Magnetization.h"
#include "OverlapTrack.h"
#include "JastrowTrack.h"
#include "StagMagnTrack.h"
#include "SlaterDeterminant.h"
#include "IdentityJastrow.h"
#include "Jastrow.h"
#include "NeelJastrowPotential.h"
#include "StagJastrowPotential.h"
#include "LogJastrowPotential.h"
#include "PairedMagnonJastrowPotential.h"
#include "RanGen.h"
#include "FileManager.h"
#include "ArgParse.h"
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <cstdlib>

using namespace std;

void myexit(int exitcode)
{
#ifdef USEMPI
    MPI_Finalize();
#endif
    exit(exitcode);
}

int main(int argc, char* argv[])
{
    // Initialize
#ifdef USEMPI
    MPI_Init(&argc,&argv);
#endif
    int comm_size(1);
    int comm_rank(0);
#ifdef USEMPI
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);
#endif
    //int wait=0;
    //if(comm_rank!=0){
    //    cout<<getpid()<<endl;
    //    wait=1;
    //    while(wait)
    //        sleep(5);
    //}
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    signal(SIGTERM,FileManager::EmergencyClose);
    Timer::tic("main");
    // calculation parameters:
    map<string,bool> bomap;
    map<string,size_t> simap;
    map<string,int> inmap;
    map<string,double> domap;
    map<string,string> stmap;
    simap["L"]=8;
    simap["samples"]=10;
    simap["samples_saves"]=1;
    simap["samples_saves_stat"]=100;
    simap["qx"]=0;
    simap["qy"]=0;
    simap["meas_interv"]=0;
    inmap["prefix"]=-1;
    inmap["therm"]=100;
    inmap["verbose"]=1;
    inmap["seed"]=time(NULL);
    domap["phi"]=0.085;
    domap["neel"]=0.055;
    domap["field"]=0.0;
    domap["neel_jastrow"]=0.0;
    domap["stag_jastrow"]=0.0;
    domap["log_jastrow"]=0.0;
    domap["pm_jastrow"]=0.0;
    domap["pm_jastrow_nn"]=0.0;
    domap["pm_jastrow_nnn"]=0.0;
    domap["phase_shift_x"]=1.0;
    domap["phase_shift_y"]=1.0;
    domap["jr"]=0.0;
    bomap["track_stagmagn"]=false;
    bomap["track_overlap"]=false;
    bomap["track_jastrow"]=false;
    bomap["meas_projheis"]=false;
    bomap["meas_stagmagn"]=false;
    bomap["meas_statspinstruct"]=false;
    bomap["meas_magnetization"]=false;
    bomap["stagflux_wav"]=false;
    bomap["sfpnxphz_wav"]=false;
    bomap["sfpnzphx_wav"]=false;
    bomap["Sztot_zero_proj"]=false;
    bomap["Sztot_non_zero_init"]=false;
    bomap["Neel_init"]=false;
    stmap["dir"]=".";
    stmap["spinstate"]="";
    stmap["channel"]="groundstate";
    int help=ArgParse(argc,argv,domap,inmap,simap,bomap,stmap);
    if(help)
        myexit(0);
    if(!comm_rank) cout<<"seed="<<inmap["seed"]<<endl;
    // Setup calculation parameters
    FileManager * fm=new FileManager(stmap["dir"],inmap["prefix"]);
    fm->Verbose()=inmap["verbose"];
    fm->MonitorTotal()=simap["samples"]*simap["samples_saves"];
    fm->StatPerSample()=simap["samples_saves"];
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
    
    RanGen::srand(inmap["seed"]+100*comm_rank);
    domap["phi"]*=M_PI;
    if(!simap["meas_interv"])
        simap["meas_interv"]=pow(simap["L"],2);

    //some input checks
    if(int(bomap["stagflux_wav"])+int(bomap["sfpnxphz_wav"])+int(bomap["sfpnzphx_wav"])!=1){
        cerr<<"Error: exactly one out of the possible trial wavefunctions [stagflux_wav,sfpnxphz_wav,sfpnzphx_wav] must be picked."<<endl;
        myexit(0);
    }

#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
    if(comm_rank){
#endif
        // Setup calculation
        double rej=0;
        size_t L=simap["L"];
        vector<size_t> Q(2);
        double neel=domap["neel"];
        double phi=domap["phi"];
        double field=domap["field"];
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
                wav=new SFpNpHxGroundState(fm,L,L,phi,neel,field,phase_shift);
            } else {//trans and long are mixed
                wav=new SFpNpHxExciton(fm,L,L,phi,neel,field,phase_shift,Q);
            }
        } else /*if(bomap["sfpnxphz_wav"])*/{
            if(stmap["channel"]=="groundstate"){
                wav=new SFpNxpHzGroundState(fm,L,L,phi,neel,field,phase_shift);
            } else {//trans and long are mixed
                wav=new SFpNxpHzExciton(fm,L,L,phi,neel,field,phase_shift,Q);
            }
        }
#ifndef NDEBUG
        cout<<"Wavefunction is:"<<endl<<*wav<<endl;
#endif
        wav->Save();
        SlaterDeterminant amp(sp,wav);
        Jastrow* jas=0;
        JastrowPotential* jaspot=0;
        if(domap["neel_jastrow"]!=0.0){
            jaspot=new NeelJastrowPotential(&slat,domap["neel_jastrow"]);
            jas=new Jastrow(sp,jaspot);
        } else if(domap["stag_jastrow"]!=0.0){
            jaspot=new StagJastrowPotential(&slat,domap["stag_jastrow"]);
            jas=new Jastrow(sp,jaspot);
        } else if(domap["log_jastrow"]!=0.0){
            jaspot=new LogJastrowPotential(&slat,domap["log_jastrow"]);
            jas=new Jastrow(sp,jaspot);
        } else if(domap["pm_jastrow"]!=0.0){
            jaspot=new PairedMagnonJastrowPotential(&slat,domap["pm_jastrow"],
                                                          domap["pm_jastrow_nn"],
                                                          domap["pm_jastrow_nnn"]);
            jas=new Jastrow(sp,jaspot);
        } else {
            jas=new IdentityJastrow;
        }
#ifndef NDEBUG
        cout<<"Initialize initial spin state"<<endl;
#endif
        if(bomap["Neel_init"]){
            vector<uint_vec_t> fst(sp->GetNfl());
            for(size_t fl=0;fl<sp->GetNfl();++fl){
                fst[fl]=uint_vec_t(sp->GetNfs()[fl],sp->GetNpt()[fl]);
            }
            for(size_t v=0;v<slat.GetNv();++v){
                bool even=(slat.GetVertices()[v]->uc[0]+slat.GetVertices()[v]->uc[1])%2;
                if(even){
                    if(bomap["stagflux_wav"])
                        fst[0][v]=0;
                    else
                        fst[0][v*2]=0;
                } else {
                    if(bomap["stagflux_wav"])
                        fst[1][v]=0;
                    else
                        fst[0][v*2+1]=0;
                }
            }
            sp->InitFock(fst);
            amp.Init();
            jas->Init();
        } else if(stmap["spinstate"]!=""){
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
                if(bomap["Sztot_non_zero_init"]){
                    int Nup=(size_t)(RanGen::uniform()*L*L);
                    pop.push_back({size_t(Nup),size_t(int(L*L)-Nup)});
                } else {
                    pop.push_back({L*L/2,L*L/2});
                }
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
                myexit(1);
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
#ifndef NDEBUG
        cout<<"initial Lattice is:"<<endl<<*sp<<endl;
        cout<<"Create Monte Carlo stepper"<<endl;
#endif
        LatticeStepper step(sp,wav,&amp,jas);
        if(!bomap["staglux_wav"] && !bomap["Sztot_zero_proj"]){
            //if applied field, Sztot is not
            //conserved and spins-up must be
            //able to tunel to spin-down.
            vector<bool> flip(1,true);
            step.SetFlavorFlip(flip);
        }
#ifndef NDEBUG
        cout<<"Create Metropolis MC object"<<endl;
#endif
        MetroMC varmc(&step,fm);
        if(bomap["meas_projheis"]){
            ProjHeis* heisen=new ProjHeis(&step,fm,&slat,domap["Bx"]);
            varmc.AddQuantity(heisen);
        }
        if(stmap["channel"]=="groundstate"){
            if(bomap["meas_stagmagn"]){
                if(bomap["stagflux_wav"]){
                    StagMagnZ* stagsz=new StagMagnZ(&step,fm);
                    varmc.AddQuantity(stagsz);
                } else {
                    StagMagn* stags=new StagMagn(&step,fm,!bomap["Sztot_zero_proj"]);
                    varmc.AddQuantity(stags);
                }
            }
            if(bomap["meas_statspinstruct"]){
                StatSpinStruct* stat=new StatSpinStruct(&step,fm,!bomap["Sztot_zero_proj"]);
                varmc.AddQuantity(stat);
            }
            if(bomap["meas_magnetization"] && (domap["Bx"]!=0 || !bomap["Sztot_conserved"])){
                Magnetization* magn=new Magnetization(&step,fm,!bomap["Sztot_zero_proj"]);
                varmc.AddQuantity(magn);
            }
        }
        if(bomap["track_stagmagn"]){
            StagMagnTrack* stt=new StagMagnTrack(&step,fm);
            varmc.AddQuantity(stt);
        }
        if(bomap["track_overlap"]){
            OverlapTrack* ovt=new OverlapTrack(&step,fm);
            varmc.AddQuantity(ovt);
        }
        if(bomap["track_jastrow"]){
            JastrowTrack* jast=new JastrowTrack(&step,fm);
            varmc.AddQuantity(jast);
        }
#ifndef NDEBUG
        cout<<"Thermalize."<<endl;
#endif
        // Start calculation: thermalize
        if(simap["therm"]){
            Timer::tic("main/thermalize");
            varmc.Walk(size_t(simap["therm"]*L*L),0);
            Timer::toc("main/thermalize");
        }
        fm->MonitorTotal()=simap["samples"]*simap["samples_saves"];         
        // Calculation
#ifndef NDEBUG
        cout<<"Random walk start now."<<endl;
#endif
        for(size_t sample=0;sample<simap["samples"];++sample){
            for(size_t m=0;m<simap["samples_saves"];++m){
                fm->MonitorCompletion()=double(sample*simap["samples_saves"]+m)/(simap["samples"]*simap["samples_saves"]);
                Timer::tic("main/ranwalk");
                varmc.Walk(simap["meas_interv"]*simap["samples_saves_stat"],simap["meas_interv"]);
                Timer::toc("main/ranwalk");
                rej=varmc.Rejection();
                for(size_t qu=0;qu<varmc.GetQuantities().size();++qu)
                    varmc.GetQuantities()[qu]->save();
                fm->DataAttribute("rej",rej);
                fm->DataAttribute("time",Timer::timer("randwalk"));
                fm->DataAttribute("statistics",(m+1)*simap["samples_saves_stat"]);
#ifdef USEMPI
                int mess(fm->message_save);
                //cout<<"rank "<<comm_rank<<": sends message_save"<<endl;
                MPI_Send(&mess,1,MPI_INT,0,fm->message_comm,MPI_COMM_WORLD);
#endif
                fm->Write();
            }
        }
#ifdef USEMPI
        int mess=fm->message_loop, stop=1;
        //cout<<"rank "<<comm_rank<<": sends message_loop"<<endl;
        MPI_Send(&mess,1,MPI_INT,0,fm->message_comm,MPI_COMM_WORLD);
        MPI_Send(&stop,1,MPI_INT,0,fm->message_loop,MPI_COMM_WORLD);
#endif
#ifndef NDEBUG
        cout<<"Calculation finished, cleaning up!"<<endl;
#endif
        sp->Save();
        delete wav;
        delete sp;
        delete jas;
        delete jaspot;
        for(size_t qu=0;qu<varmc.GetQuantities().size();++qu)
            delete varmc.GetQuantities()[qu];
#ifdef USEMPI
    } else {
        fm->MainLoop();
    }
#endif
    // Output
    Timer::toc("main");
    if(comm_rank==0){
        for(int r=0;r<comm_size;++r){
            std::cout<<"############################################"<<std::endl;
            std::cout<<"rank="<<r<<endl;
            if(r==0)
                std::cout<<Timer::report();
#ifdef USEMPI
            else {
                int len(0);
                MPI_Recv(&len,1,MPI_INT,r,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                char *mes=new char[len+1];
                MPI_Recv(mes,len+1,MPI_CHAR,r,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                cout<<mes<<endl;
                delete [] mes;
            }
#endif
            std::cout<<"############################################"<<std::endl;
        }
        std::cout<<"output prefix="<<fm->Prefix()<<std::endl;
    } else {
#ifdef USEMPI
        string rep=Timer::report();
        int len=rep.size();
        MPI_Send(&len,1,MPI_INT,0,1,MPI_COMM_WORLD);
        char* mes=new char[len+1];
        memcpy(mes,rep.c_str(),(len+1)*sizeof(char));
        MPI_Send(mes,len+1,MPI_CHAR,0,1,MPI_COMM_WORLD);
        delete [] mes;
#endif
    }
    delete fm;
#ifdef USEMPI
    MPI_Finalize();
#endif
    return 0;
}

